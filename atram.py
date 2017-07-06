"""The aTRAM assembly program."""

import os
import re
import sys
import math
import argparse
import textwrap
import tempfile
import datetime
import subprocess
import multiprocessing
from shutil import which
import psutil
from Bio import SeqIO
import lib.db as db
import lib.bio as bio
import lib.log as log
import lib.blast as blast
import lib.file_util as file_util
from lib.assembler import factory


def run(args):
    """Setup and run atram."""

    log.setup(args)

    all_shards = fraction_of_shards(args)

    db_conn = db.connect(args.blast_db)

    # Make sure the database version matches what we built it with
    version = db.get_version(db_conn)
    if version != db.VERSION:
        log.fatal('The database was built with version {} but you are running '
                  'version {}. You need to rebuild the atram database by '
                  'running atram_preprocessor.py again.'.format(
                      version, db.VERSION))

    assembler = factory(args) if args.assembler else None
    query = initialize_query(args, db_conn, assembler)

    atram_loop(args, db_conn, assembler, query, all_shards)

    if assembler:
        output_assembly_results(args, db_conn)
    else:
        output_blast_only_results(args, db_conn)

    db_conn.close()


def atram_loop(args, db_conn, assembler, query, all_shards):
    """The main program loop."""

    for iteration in range(args.start_iteration, args.iterations + 1):
        log.info('aTRAM iteration %i' % iteration, line_break='')

        blast_query_against_all_sras(args, query, all_shards, iteration)

        # If we don't have an assembler then we just want the blast hits
        if not args.assembler:
            break

        # Exit if there are no blast hits
        if not db.sra_blast_hits_count(db_conn, iteration):
            log.info('No blast hits in iteration %i' % iteration)
            break

        assembler.initialize_iteration(iteration)

        assembler.write_input_files(db_conn)

        try:
            log.info('Assembling shards with {}: iteration {}'.format(
                args.assembler, iteration))
            assembler.assemble()
        except TimeoutError:
            msg = 'Time ran out for the assembler after {} (HH:MM:SS)'.format(
                datetime.timedelta(seconds=args.timeout))
            log.fatal(msg)
        except subprocess.CalledProcessError as cpe:
            msg = 'The assembler failed with error: ' + str(cpe)
            log.fatal(msg)

        # Exit if nothing was assembled
        if not os.path.exists(assembler.file['output']) \
                or not os.path.getsize(assembler.file['output']):
            log.info('No new assemblies in iteration %i' % iteration)
            break

        high_score = filter_contigs(args, db_conn, assembler, iteration)

        count = db.assembled_contigs_count(
            db_conn, iteration, args.bit_score, args.contig_length)
        if not count:
            log.info(('No contigs had a bit score greater than {} and are at '
                      'least {} long in iteration {}. The highest score for '
                      'this iteration is {}').format(args.bit_score,
                                                     args.contig_length,
                                                     iteration,
                                                     high_score))
            break

        if count == db.iteration_overlap_count(
                db_conn, iteration, args.bit_score, args.contig_length):
            log.info('No new contigs were found in iteration %i' % iteration)
            break

        query = create_queries_from_contigs(
            db_conn, assembler, iteration, args.bit_score, args.contig_length)
    else:
        log.info('All iterations completed', line_break='')


def initialize_query(args, db_conn, assembler):
    """Get the first set of query sequences."""

    if args.start_iteration < 2:
        db.create_sra_blast_hits_table(db_conn)
        db.create_contig_blast_hits_table(db_conn)
        db.create_assembled_contigs_table(db_conn)
        query = args.query
    else:
        query = create_queries_from_contigs(db_conn,
                                            assembler,
                                            args.start_iteration - 1,
                                            args.bit_score,
                                            args.contig_length)

    if os.path.getsize(query) < 10:
        err = 'There are no sequences '
        if args.start_iteration < 2:
            err += 'in {}'.format(args.query)
        else:
            err += 'for starting iteration {}'.format(args.start_iteration)
        log.fatal(err)

    return query


def fraction_of_shards(args):
    """Get the fraction of shards to use."""

    all_shards = blast.all_shard_paths(args.blast_db)
    last_index = int(len(all_shards) * args.fraction)

    return all_shards[:last_index]


def blast_query_against_all_sras(args, query, all_shards, iteration):
    """Blast the queries against the SRA databases. We're using a
    map-reduce strategy here. We map the blasting of the query sequences
    and reduce the output into one fasta file.
    """

    log.info('Blasting query against shards: iteration %i' % iteration)

    # for shard_path in all_shards:
    #     blast_query_against_sra(
    #         dict(vars(args)), shard_path, query, iteration)
    with multiprocessing.Pool(processes=args.cpus) as pool:
        results = [pool.apply_async(
            blast_query_against_sra,
            (dict(vars(args)), shard_path, query, iteration))
                   for shard_path in all_shards]
        _ = [result.get() for result in results]


def blast_query_against_sra(args, shard_path, query, iteration):
    """Blast the query against one blast DB shard. Then write the results to
    the database.
    """
    # NOTE: Because this is called in a child process, the address space is not
    # shared with the parent (caller) hence we cannot share object variables.

    output_file = blast.output_file_name(
        args['temp_dir'], shard_path, iteration)

    blast.against_sra(args, shard_path, query, output_file, iteration)

    db_conn = db.connect(args['blast_db'])

    shard = os.path.basename(shard_path)

    batch = []

    for hit in blast.hits(output_file):
        match = blast.PARSE_RESULTS.match(hit['title'])
        if match:
            seq_name = match.group(1)
            seq_end = match.group(2)
        else:
            seq_name = hit['title']
            seq_end = ''
        batch.append((iteration, seq_end, seq_name, shard))
    db.insert_blast_hit_batch(db_conn, batch)
    db_conn.close()


def output_blast_only_results(args, db_conn):
    """Output this file if we are not assembling the contigs."""

    log.info('Output blast only results')

    file_name = file_util.output_file(args, 'blast_only.fasta')
    with open(file_name, 'w') as out_file:
        for row in db.get_sra_blast_hits(db_conn, 1):
            out_file.write('>{}{}\n'.format(row['seq_name'], row['seq_end']))
            out_file.write('{}\n'.format(row['seq']))


def filter_contigs(args, db_conn, assembler, iteration):
    """Remove junk from the assembled contigs."""

    log.info('Saving assembled contigs: iteration %i' % iteration)

    blast_db = blast.temp_db_name(args.temp_dir, args.blast_db, iteration)
    hits_file = blast.output_file_name(args.temp_dir, args.blast_db, iteration)

    blast.create_db(args.temp_dir, assembler.file['output'], blast_db)

    blast.against_contigs(args, blast_db, args.query, hits_file)

    save_blast_against_contigs(db_conn, assembler, hits_file, iteration)

    all_hits = {row['contig_id']: row
                for row
                in db.get_contig_blast_hits(db_conn, iteration)}

    return save_contigs(db_conn, all_hits, assembler, iteration)


def save_blast_against_contigs(db_conn, assembler, hits_file, iteration):
    """Only save contigs that have bit scores above the cut-off."""

    batch = []

    for hit in blast.hits(hits_file):
        contig_id = assembler.parse_contig_id(hit['title'])
        batch.append((iteration, contig_id, hit['title'],
                      hit['bit_score'], hit['len'],
                      hit['query_from'], hit['query_to'],
                      hit.get('query_strand', ''),
                      hit['hit_from'], hit['hit_to'],
                      hit.get('hit_strand', '')))

    db.insert_contig_hit_batch(db_conn, batch)


def save_contigs(db_conn, all_hits, assembler, iteration):
    """Save the contigs to the database."""

    batch = []
    high_score = 0
    with open(assembler.file['output']) as in_file:
        for contig in SeqIO.parse(in_file, 'fasta'):
            contig_id = assembler.parse_contig_id(contig.description)
            if contig_id in all_hits:
                hit = all_hits[contig_id]
                batch.append((
                    iteration, contig.id, str(contig.seq), contig.description,
                    hit['bit_score'], hit['len'],
                    hit['query_from'], hit['query_to'], hit['query_strand'],
                    hit['hit_from'], hit['hit_to'], hit['hit_strand']))
    db.insert_assembled_contigs_batch(db_conn, batch)

    return high_score


def create_queries_from_contigs(
        db_conn, assembler, iteration, bit_score, length):
    """Crate a new file with the contigs that will be used as the
    next query query.
    """

    log.info('Creating new query files: iteration %i' % iteration)

    query = assembler.path('long_reads.fasta')
    assembler.file['long_reads'] = query

    with open(query, 'w') as query_file:
        for row in db.get_assembled_contigs(
                db_conn, iteration, bit_score, length):
            query_file.write('>{}\n'.format(row[0]))
            query_file.write('{}\n'.format(row[1]))

    return query


def output_assembly_results(args, db_conn):
    """Write the assembled contigs to a fasta file."""

    file_name = file_util.output_file(args, 'filtered_contigs.fasta')
    with open(file_name, 'w') as out_file:
        for row in db.get_all_assembled_contigs(
                db_conn, args.bit_score, args.contig_length):
            output_one_assembly(out_file, row)

    file_name = file_util.output_file(args, 'all_contigs.fasta')
    with open(file_name, 'w') as out_file:
        for row in db.get_all_assembled_contigs(db_conn):
            output_one_assembly(out_file, row)


def output_one_assembly(out_file, row):
    """Write an assembly to the output fasta file."""

    seq = row['seq']
    suffix = ''
    if row['query_strand'] and row['hit_strand'] and \
            row['query_strand'] != row['hit_strand']:
        seq = bio.reverse_complement(seq)
        suffix = '_REV'

    header = ('>{}_{}{} iteration={} contig_id={} '
              'score={}\n').format(
                  row['iteration'], row['contig_id'], suffix,
                  row['iteration'], row['contig_id'], row['bit_score'])
    out_file.write(header)
    out_file.write('{}\n'.format(seq))


def parse_command_line(temp_dir):
    """Process command-line arguments."""

    description = """
        This is the aTRAM script. It takes a query sequence and
        a set of blast databases built with the atram_preprocessor.py script
        and builds an assembly.
        """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(description))

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(db.VERSION))

    required_args(parser)
    optional_atram_args(parser)
    optional_blast_filter_args(parser)
    optional_blast_args(parser)
    optional_assembler_args(parser)

    args = parser.parse_args()

    check_required_args(args)
    check_optional_atram_args(args)
    check_optional_blast_args(args, temp_dir)
    check_optional_assembler_args(args)

    find_programs(args)

    return args


def required_args(parser):
    """Add required arguments to the parser."""

    group = parser.add_argument_group('required arguments')

    group.add_argument('-b', '--blast-db', '--sra', '--db', '--database',
                       required=True, metavar='DB',
                       help='This needs to match the DB prefix you '
                            'entered for atram_preprocessor.py.')

    group.add_argument('-o', '--output', required=True,
                       help='This is the prefix of all of the output files. '
                            'So you can identify different blast output file '
                            'sets. You may include a directory as part of the '
                            'prefix.')

    group.add_argument('-q', '--query', '--target',
                       help='The path to the fasta file with sequences of '
                            'interest. Required unless you specify a '
                            '"--start-iteration".')


def check_required_args(args):
    """Make sure required arguments are reasonable."""

    # Touch up blast DB name
    pattern = (r'^ (.*?)'
               r'(  \.atram(_preprocessor)?\.log'
               r' | \.blast_\d{3}\.(nhr|nin|nsq)'
               r' | \.sqlite\.db  )?$')
    args.blast_db = re.sub(pattern, r'\1', args.blast_db, re.I | re.X)

    # Check query
    if not args.query and not args.start_iteration:
        log.fatal('We need a "--query" sequence.')


def optional_atram_args(parser):
    """Add optional atram arguments to the parser."""

    group = parser.add_argument_group('optional aTRAM arguments')

    group.add_argument('-a', '--assembler',
                       choices=['abyss', 'trinity', 'velvet', 'spades'],
                       help='Which assembler to use. If you do not use this '
                            'argument then aTRAM will do a single blast run '
                            'and stop before assembly.')

    group.add_argument('-i', '--iterations', type=int, default=5, metavar='N',
                       help='The number of pipeline iterations. '
                            'The default is "5".')

    group.add_argument('-p', '--protein', action='store_true',
                       help='Are the query sequences protein? '
                            'aTRAM will guess if you skip this argument.')

    group.add_argument('--fraction', type=float, default=1.0,
                       help='Use only the specified fraction of the aTRAM '
                            'database. The default is "1.0"')

    # group.add_argument('--complete', action='store_true',
    #                    help='Automatically quit when a complete homolog is '
    #                         'recovered.')

    cpus = min(10, os.cpu_count() - 4 if os.cpu_count() > 4 else 1)
    group.add_argument('--cpus', '--processes', '--max-processes',
                       type=int, default=cpus,
                       help=('Number of cpus to use. This will also be used '
                             'for the assemblers when possible. Defaults to: '
                             'Total CPUs - 4 = "{}"').format(cpus))

    group.add_argument('--log-file',
                       help='Log file (full path). The default is to use the '
                            'DIR and DB arguments to come up with a name like '
                            'so "DIR/DB_atram.log"')

    group.add_argument('--path',
                       help='If the assembler or blast you want to use is not '
                            'in your $PATH then use this to prepend '
                            'directories to your path.')

    group.add_argument('--start-iteration', '--restart',
                       type=int, default=1, metavar='N',
                       help='If resuming from a previous run, which iteration '
                            'number to start from. The default is "1".')

    group.add_argument('-t', '--temp-dir', metavar='DIR',
                       help='You may save intermediate files for debugging '
                            'in this directory. The directory must be empty.')

    group.add_argument('-T', '--timeout', metavar='SECONDS', default=300,
                       type=int,
                       help='How many seconds to wait for an assembler before '
                            'stopping the run. To wait forever set this to 0. '
                            'The default is "300" (5 minutes).')


def check_optional_atram_args(args):
    """Make sure optional atram arguments are reasonable."""

    # Set default log file name
    if not args.log_file:
        args.log_file = '{}.{}.log'.format(
            args.blast_db, os.path.basename(sys.argv[0][:-3]))

    # If not --protein then probe to see if it's a protein seq
    if not args.protein:
        with open(args.query) as in_file:
            for query in SeqIO.parse(in_file, 'fasta'):
                if bio.is_protein(str(query.seq)):
                    args.protein = True

    # Prepend to PATH environment variable if requested
    if args.path:
        os.environ['PATH'] = '{}:{}'.format(args.path, os.environ['PATH'])


def optional_blast_filter_args(parser):
    """Add optional values for blast-filtering contigs to the parser."""

    group = parser.add_argument_group(
        'optional values for blast-filtering contigs')

    group.add_argument('--bit-score', type=float, default=70.0,
                       metavar='SCORE',
                       help='Remove contigs that have a value less than this. '
                            'The default is "70.0"')

    group.add_argument('--contig-length', '--length', type=int, default=100,
                       help='Remove blast hits that are shorter than this '
                            'length. The default is "100".')


def optional_blast_args(parser):
    """Add optional blast arguments to the parser."""

    group = parser.add_argument_group('optional blast arguments')

    group.add_argument('--db-gencode', type=int, default=1,
                       metavar='CODE',
                       help='The genetic code to use during blast runs. '
                            'The default is "1".')

    group.add_argument('--evalue', type=float, default=1e-10,
                       help='The default evalue is "1e-10".')

    group.add_argument('--max-target-seqs', type=int, default=100000000,
                       metavar='MAX',
                       help='Maximum hit sequences per shard. '
                            'Default is calculated based on the available '
                            'memory and the number of shards. ')


def check_optional_blast_args(args, temp_dir):
    """Make sure optional blast arguments are reasonable."""

    # Calculate the default max_target_seqs per shard
    if not args.max_target_seqs:
        all_shards = blast.all_shard_paths(args.blast_db)
        args.max_target_seqs = int(2 * args.max_memory / len(all_shards)) * 1e6

    file_util.temp_root_dir(args, temp_dir)


def optional_assembler_args(parser):
    """Add optional assembler arguments to the parser."""

    group = parser.add_argument_group('optional assembler arguments')

    group.add_argument('--no-long-reads', action='store_true',
                       help='Do not use long reads during assembly. '
                            '(Abyss, Trinity, Velvet)')

    group.add_argument('--kmer', type=int, default=64,
                       help='k-mer size. The default is "64" for Abyss and '
                            '"31" for Velvet. Note: the maximum kmer length '
                            'for Velvet is 31. (Abyss, Velvet)')

    group.add_argument('--mpi', action='store_true',
                       help='Use MPI for this assembler. The assembler '
                            'must have been compiled to use MPI. (Abyss)')

    group.add_argument('--bowtie2', action='store_true',
                       help='Use bowtie2 during assembly. (Trinity)')

    max_mem = max(1.0, math.floor(
        psutil.virtual_memory().available / 1024**3 / 2))
    group.add_argument('--max-memory',
                       default=max_mem, metavar='MEMORY', type=int,
                       help='Maximum amount of memory to use in gigabytes. '
                            'The default is "{}". (Trinity, Spades)'.format(
                                max_mem))

    group.add_argument('--exp-coverage', '--expected-coverage',
                       type=int, default=30,
                       help='The expected coverage of the region. '
                            'The default is "30". (Velvet)')

    group.add_argument('--ins-length', type=int, default=300,
                       help='The size of the fragments used in the short-read '
                            'library. The default is "300". (Velvet)')

    group.add_argument('--min-contig-length', type=int, default=100,
                       help='The minimum contig length used by the assembler '
                            'itself. The default is "100". (Velvet)')

    group.add_argument('--cov-cutoff', default='off',
                       help='Read coverage cutoff value. Must be a positive '
                            'float value, or "auto", or "off". '
                            'The default is "off". (Spades)')


def check_optional_assembler_args(args):
    """Make sure optional assembler arguments are reasonable."""

    # Check kmer
    if args.assembler == 'velvet' and args.kmer > 31:
        args.kmer = 31

    # Check cov_cutoff
    if args.cov_cutoff not in ['off', 'auto']:
        err = ('Read coverage cutoff value. Must be a positive '
               'float value, or "auto", or "off"')
        try:
            value = float(args.cov_cutoff)
        except ValueError:
            log.fatal(err)
        if value < 0:
            log.fatal(err)


def find_programs(args):
    """Make sure we can find the programs needed by the assembler and blast."""

    if not (which('makeblastdb') and which('tblastn') and which('blastn')):
        err = ('We could not find the programs "makeblastdb", "tblastn", or '
               '"blastn". You either need to install them or you need adjust '
               'the PATH environment variable with the "--path" option so '
               'that aTRAM can find it.')
        log.fatal(err)

    if args.assembler == 'abyss' and not which('abyss-pe'):
        err = ('We could not find the "abyss-pe" program. You either need to '
               'install it or you need to adjust the PATH environment '
               'variable with the "--path" option so that aTRAM can find it.')
        log.fatal(err)

    if args.assembler == 'abyss' and not args.no_long_reads \
            and not which('bwa'):
        err = ('We could not find the "bwa-mem" program. You either need to '
               'install it, adjust the PATH environment variable '
               'with the "--path" option, or you may use the '
               '"--no-long-reads" option to not use this program.')
        log.fatal(err)

    if args.assembler == 'trinity' and not which('Trinity'):
        err = ('We could not find the "Trinity" program. You either need to '
               'install it or you need to adjust the PATH environment '
               'variable with the "--path" option so that aTRAM can find it.')
        log.fatal(err)

    if args.assembler == 'trinity' and args.bowtie2 and not which('bowtie2'):
        err = ('We could not find the "bowtie2" program. You either need to '
               'install it, adjust the PATH environment variable '
               'with the "--path" option, or you may skip using this program '
               'by not using the "--bowtie2" option.')
        log.fatal(err)

    if args.assembler == 'velvet' and \
            not (which('velveth') and which('velvetg')):
        err = ('We could not find either the "velveth" or "velvetg" program. '
               'You either need to install it or you need to adjust the PATH '
               'environment variable with the "--path" option so that aTRAM '
               'can find it.')
        log.fatal(err)

    if args.assembler == 'spades' and not which('spades.py'):
        err = ('We could not find the "Spades" program. You either need to '
               'install it or you need to adjust the PATH environment '
               'variable with the "--path" option so that aTRAM can find it.')
        log.fatal(err)


if __name__ == '__main__':

    with tempfile.TemporaryDirectory(prefix='atram_') as TEMP_DIR:
        ARGS = parse_command_line(TEMP_DIR)
        run(ARGS)
