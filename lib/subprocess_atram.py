"""Atram subprocess functions."""

from os.path import basename
import lib.db as db
import lib.db_atram as db_atram
import lib.blast as blast


def blast_query(args, state, shard):
    """
    Blast the query against one blast DB shard.

    Then write the results to a temporary database.
    """
    output_file = blast.output_file_name(state['iter_dir'], shard)

    blast.against_sra(args, state, output_file, shard)

    shard = basename(shard)

    batch = []

    with db.connect(state['blast_db']) as cxn:
        db.aux_db(
            cxn,
            args['temp_dir'],
            state['blast_db'],
            state['query_target'])
        is_single_end = db.is_single_end(cxn)
        db.aux_detach(cxn)

    hits = blast.hits(output_file)
    for hit in hits:
        seq_name, seq_end = blast.parse_blast_title(
            hit['title'], is_single_end)
        batch.append((state['iteration'], seq_end, seq_name, shard))

    with db.subprocess_db(
            args['temp_dir'],
            state['blast_db'],
            state['query_target'],
            shard) as cxn:
        db_atram.insert_blast_hit_batch(cxn, batch)
