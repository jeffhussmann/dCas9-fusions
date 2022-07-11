import argparse
import logging
import sys
from pathlib import Path

import tqdm

import knock_knock.experiment
import knock_knock.parallel

import dCas9_fusions.experiment

def parallel(args):
    logger = logging.getLogger(__name__)
    logger.propagate = False
    logger.setLevel(logging.DEBUG)
    stream_handler = logging.StreamHandler()
    formatter = logging.Formatter(fmt='%(asctime)s: %(message)s',
                                  datefmt='%y-%m-%d %H:%M:%S',
                                 )
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.DEBUG)
    logger.addHandler(stream_handler)

    if args.group:
        args.conditions['batch'] = args.group

    exps = knock_knock.experiment.get_all_experiments(args.project_directory, args.conditions)

    def process_stage(stage):
        with knock_knock.parallel.PoolWithLoggerThread(args.max_procs, logger) as process_pool:
            arg_tuples = []

            for _, exp in exps.items():
                arg_tuple = (exp.base_dir, exp.batch, exp.sample_name, stage, args.progress, True)
                arg_tuples.append(arg_tuple)

            process_pool.starmap(knock_knock.experiment.process_experiment_stage, arg_tuples)

    stages = args.stages.split(',')
    for stage in stages:
        process_stage(stage)

def process(args):
    stages = args.stages.split(',')

    sample_sheet = knock_knock.experiment.load_sample_sheet(args.base_dir, args.batch_name)
    description = sample_sheet[args.sample_name]

    exp = dCas9_fusions.experiment.dCas9FusionExperiment(args.base_dir,
                                                         args.batch_name,
                                                         args.sample_name,
                                                         description=description,
                                                        )

    for stage in stages:
        logging.info(f'Starting {exp} {stage}')
        exp.process(stage)

if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s: %(message)s',
                        datefmt='%y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                       )

    parser = argparse.ArgumentParser(prog='count_domains')

    subparsers = parser.add_subparsers(dest='subcommand', title='subcommands')
    subparsers.required = True

    def add_project_directory_arg(parser):
        parser.add_argument('base_dir', type=Path, help='the base directory to store input data, reference annotations, and analysis output for a project')

    parser_process = subparsers.add_parser('process', help='process a single sample')
    add_project_directory_arg(parser_process)
    parser_process.add_argument('batch_name', help='batch name')
    parser_process.add_argument('sample_name', help='sample name')
    parser_process.add_argument('--progress', const=tqdm.tqdm, action='store_const', help='show progress bars')
    parser_process.add_argument('--stages', default='preprocess,align,categorize,visualize')
    parser_process.set_defaults(func=process)

    parser_parallel = subparsers.add_parser('parallel', help='process multiple samples in parallel')
    add_project_directory_arg(parser_parallel)
    parser_parallel.add_argument('max_procs', type=int, help='maximum number of samples to process at once')
    parser_parallel.add_argument('--batch', help='if specified, the single batch name to process; if not specified, all groups will be processed')
    parser_parallel.add_argument('--stages', default='preprocess,align,categorize,visualize')
    parser_parallel.add_argument('--progress', const=tqdm.tqdm, action='store_const', help='show progress bars')
    parser_parallel.set_defaults(func=parallel)

    args = parser.parse_args()
    args.func(args)