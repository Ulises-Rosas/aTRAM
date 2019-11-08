"""Build assemblies using the aTRAM algorithm.."""

import re
import os
from os.path import basename, split, splitext, join
from multiprocessing import Pool
from Bio import SeqIO
import lib.db as db
import lib.db_atram as db_atram
import lib.subprocess_atram as subprocess_atram
import lib.log as log
import lib.bio as bio
import lib.util as util
import lib.blast as blast
import lib.assembler as assembler


def assemble(args):
    """Loop thru every blast/query pair and run an assembler for each one."""
    with util.make_temp_dir(
            where=args['temp_dir'],
            prefix='atram_',
            keep=args['keep_temp_dir']) as temp_dir:
        util.update_temp_dir(temp_dir, args)

        queries = split_queries(args)

        for blast_db in args['blast_db']:

            with db.connect(blast_db, check_version=True) as cxn:
                for query in queries:
                    db.aux_db(cxn, args['temp_dir'], blast_db, query)
                    clean_database(cxn)

                    log.setup(args['log_file'], blast_db, query)

                    asm = assembler.factory(args, cxn)

                    try:
                        assembly_loop(args, asm, blast_db, query)
                    except (TimeoutError, RuntimeError):
                        pass
                    except Exception as err:  # pylint: disable=broad-except
                        log.error('Exception: {}'.format(err))
                    finally:
                        asm.write_final_output(blast_db, query)

                    db.aux_detach(cxn)


def assembly_loop(args, asm, blast_db, query):
    """Iterate over the assembly processes."""
    for iteration in range(1, asm.args['iterations'] + 1):
        log.info('aTRAM blast DB = "{}", query = "{}", iteration {}'.format(
            blast_db, split(query)[1], iteration))

        asm.init_iteration(blast_db, query, iteration)

        with util.make_temp_dir(
                where=args['temp_dir'],
                prefix=asm.file_prefix(),
                keep=args['keep_temp_dir']) as iter_dir:

            asm.setup_files(iter_dir)

            query = assembly_loop_iteration(args, asm)

            if not query:
                break

    else:
        log.info('All iterations completed')


def assembly_loop_iteration(args, asm):
    """One iteration of the assembly loop."""
    blast_query(asm)

    count = asm.count_blast_hits()
    if asm.blast_only or count == 0:
        return False

    asm.write_input_files()

    asm.run()

    if asm.nothing_assembled():
        return False

    high_score = filter_contigs(asm)

    count = asm.assembled_contigs_count(high_score)

    if not count:
        return False

    if asm.no_new_contigs(count):
        return False

    return create_query_from_contigs(args, asm)


def split_queries(args):
    """
    Create query target for every query and query-split file.

    We put each query record into its own file for blast queries.
    """
    if not args.get('query_split'):
        return args['query'][:]

    queries = []

    path = join(args['temp_dir'], 'queries')
    os.makedirs(path, exist_ok=True)

    for query_path in args['query_split']:

        query_name = splitext(basename(query_path))[0]

        with open(query_path) as query_file:
            for i, rec in enumerate(SeqIO.parse(query_file, 'fasta'), 1):

                query_id = re.sub(r'\W+', '_', rec.id)

                query_file = join(
                    path,
                    '{}_{}_{}.fasta'.format(query_name, query_id, i))

                write_query_seq(query_file, rec.id, str(rec.seq))

                queries.append(query_file)

    if not args.get('protein'):
        args['protein'] = bio.fasta_file_has_protein(queries)

    return queries


def write_query_seq(file_name, seq_id, seq):
    """Write the sequence to a fasta file."""
    with open(file_name, 'w') as query_file:
        util.write_fasta_record(query_file, seq_id, seq)


def clean_database(cxn):
    """Create database tables for an atram run."""
    db_atram.create_sra_blast_hits_table(cxn)
    db_atram.create_contig_blast_hits_table(cxn)
    db_atram.create_assembled_contigs_table(cxn)


def clean_subprocess_databases(asm, all_shards):
    """Create database tables for the subprocess databases."""
    if asm.state['iteration'] > 1:
        return

    for shard in all_shards:
        with db.subprocess_db(
                asm.args['temp_dir'],
                asm.state['blast_db'],
                asm.state['query_target'],
                shard) as cxn:
            db_atram.create_sra_blast_hits_table(cxn, 'subprocess')


def blast_query(asm):
    """
    Blast the query against the SRA databases.

    We're using a map-reduce strategy here. We map the blasting of the query
    sequences and reduce the output into one fasta file.
    """
    log.info('Blasting query against shards: iteration {}'.format(
        asm.state['iteration']))

    all_shards = shard_fraction(asm)
    clean_subprocess_databases(asm, all_shards)

    with Pool(processes=asm.args['cpus']) as pool:
        results = [pool.apply_async(
            subprocess_atram.blast_query,
            (asm.args, asm.simple_state(), shard))
                   for shard in all_shards]
        all_results = [result.get() for result in results]

    transfer_blast_hits(asm, all_shards)

    log.info('All {} blast results completed'.format(len(all_results)))


def transfer_blast_hits(asm, all_shards):
    """Copy blast results from shard databases into the main auxiliary DB."""
    cxn = asm.state['cxn']
    for shard in all_shards:
        db.subprocess_attach(
            cxn,
            asm.args['temp_dir'],
            asm.state['blast_db'],
            asm.state['query_target'],
            shard)
        db_atram.transfer_blast_hits(cxn, asm.state['iteration'])
        db.subprocess_detach(cxn)


def shard_fraction(asm):
    """
    Get the shards we are using.

    We may not want the entire DB for highly redundant libraries.
    """
    all_shards = blast.all_shard_paths(asm.state['blast_db'])
    last_index = int(len(all_shards) * asm.args['fraction'])
    return all_shards[:last_index]


def filter_contigs(asm):
    """Remove junk from the assembled contigs."""
    log.info('Saving assembled contigs: iteration {}'.format(
        asm.state['iteration']))

    blast_db = blast.temp_db_name(
        asm.state['iter_dir'], asm.state['blast_db'])

    hits_file = blast.output_file_name(
        asm.state['iter_dir'], asm.state['blast_db'])

    blast.create_db(
        asm.state['iter_dir'], asm.file['output'], blast_db)

    blast.against_contigs(
        blast_db,
        asm.state['query_target'],
        hits_file,
        protein=asm.args['protein'],
        db_gencode=asm.args['db_gencode'],
        temp_dir=asm.args['temp_dir'],
        timeout=asm.args['timeout'])

    save_blast_against_contigs(asm, hits_file)

    all_hits = {row['contig_id']: row
                for row
                in db_atram.get_contig_blast_hits(
                    asm.state['cxn'],
                    asm.state['iteration'])}

    return save_contigs(asm, all_hits)


def save_blast_against_contigs(asm, hits_file):
    """Save all of the blast hits."""
    batch = []

    for hit in blast.hits(hits_file):
        contig_id = asm.parse_contig_id(hit['title'])
        batch.append((
            asm.state['iteration'],
            contig_id,
            hit['title'],
            hit['bit_score'],
            hit['len'],
            hit['query_from'],
            hit['query_to'],
            hit.get('query_strand', ''),
            hit['hit_from'],
            hit['hit_to'],
            hit.get('hit_strand', '')))

    db_atram.insert_contig_hit_batch(asm.state['cxn'], batch)


def save_contigs(asm, all_hits):
    """Save the contigs to the database."""
    batch = []
    high_score = 0
    with open(asm.file['output']) as in_file:
        for contig in SeqIO.parse(in_file, 'fasta'):
            contig_id = asm.parse_contig_id(contig.description)
            if contig_id in all_hits:
                hit = all_hits[contig_id]
                batch.append((
                    asm.state['iteration'],
                    contig.id,
                    str(contig.seq),
                    contig.description,
                    hit['bit_score'],
                    hit['len'],
                    hit['query_from'],
                    hit['query_to'],
                    hit['query_strand'],
                    hit['hit_from'],
                    hit['hit_to'],
                    hit['hit_strand']))
    db_atram.insert_assembled_contigs_batch(asm.state['cxn'], batch)

    return high_score


def create_query_from_contigs(args, asm):
    """Crate a new file with the contigs used as the next query."""
    log.info('Creating new query files: iteration {}'.format(
        asm.state['iteration']))

    query_dir = join(args['temp_dir'], 'queries')
    os.makedirs(query_dir, exist_ok=True)

    query_file = asm.file_prefix() + 'long_reads.fasta'
    query = join(query_dir, query_file)
    asm.file['long_reads'] = query

    with open(query, 'w') as query_file:
        for row in db_atram.get_assembled_contigs(
                asm.state['cxn'],
                asm.state['iteration'],
                asm.args['bit_score'],
                asm.args['contig_length']):
            util.write_fasta_record(query_file, row[0], row[1])

    return query
