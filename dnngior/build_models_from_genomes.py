#!/usr/bin/env python3

__author__ = 'M.D. Boer'
__version__ = '0.1'
__date__ = '15 Jan, 2024'

import argparse
import os
import sys

def parse_arguments():
    class PathAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            path = os.path.expanduser(values.rstrip('/'))

            if not path.startswith('/') and not path.startswith('.'):
                path = f'./{path}'

            setattr(namespace, self.dest, path)

    parser = argparse.ArgumentParser(
            prog='build_models_from_genome.py',
            description=('Command line script to build models from a folder with fasta files (-f DIR) or '
            'a folder with ungapfilled base models (-m DIR)'),
            usage=('python build_model_from_genome.py (-f [DIR] | -d [DIR]) -o output_folder '),
            add_help=False
            )

    required_choice = parser.add_argument_group('Required choice')
    group = required_choice.add_mutually_exclusive_group(required=True)
    group.add_argument(
            '-f',
            '--fasta_files',
            dest='fasta_folder',
            metavar='',
            type=str,
            action=PathAction,
            help=(
                '[DIR] '
                'folder containing the fasta files you want to build models for'
                )
            )
    group.add_argument(
            '-m',
            '--models',
            dest='model_folder',
            metavar='',
            type=str,
            action=PathAction,
            help=(
                '[DIR] '
                'Folder in case you allready have base models'
            )
    )

    required = parser.add_argument_group('Required')
    required.add_argument(
            '-o',
            '--output',
            dest='output_folder',
            metavar='',
            required=True,
            type=str,
            action=PathAction,
            help=('[DIR] Path to output folder.')
            )

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument(
            '-v',
            '--version',
            action='version',
            version=(f'v{__version__} ({__date__}).'),
            help='Print version information and exit.'
            )
    optional.add_argument(
            '-h',
            '--help',
            action='help',
            help='Show this help message and exit.'
            )

    (args, extra_args) = parser.parse_known_args()
    if len(extra_args) > 1:
        sys.exit('error: to many arguments supplied:\n{0}'.format(
            '\n'.join(extra_args)))

    return args

def lazy_loading_modules():
    from dnngior import Gapfill, NN_Predictor
    from dnngior.reaction_class import Reaction
    from modelseedpy import MSBuilder, MSGenome
    from modelseedpy.core import msmedia
    from modelseedpy.core.rast_client import RastClient
    from cobra.io import write_sbml_model

def build_base_model(path_to_fasta):
    # Set the path to your genome
    print("Build MSGenome object")
    patric_genome = MSGenome.from_fasta(path_to_fasta, split = ' ')

    rast = RastClient()
    rast.annotate_genome(patric_genome)

    print("Start building base model")
    base_model = MSBuilder.build_metabolic_model(model_id = os.path.basename(path_to_fasta),
                                                 genome   = patric_genome,
                                                 index    = "0",
                                                 classic_biomass = True,
                                                 gapfill_model   = False,
                                                 gapfill_media   = None,
                                                 annotate_with_rast = False,
                                                 allow_all_non_grp_reactions = True
                                            )
    print('base model created')
    return base_model

def build_gapfilled_model_from_fasta(path_to_fasta, args):
    model_name = os.path.basename(path_to_fasta)[:-4]
    ug_model = build_base_model(path_to_fasta)
    print('#reactions: {}'.format(len(ug_model.reactions.list_attr('id'))))
    ug_location = os.path.join(args.output_folder,'base_models','base_{}'.format(model_name))
    write_sbml_model(cobra_model = ug_model, filename =ug_location)
    gapfill_model_wrapper(ug_location, args)

def gapfill_model_wrapper(ug_location, args):
    gf_model = Gapfill(ug_location).gapfilledModel
    ug_model_name = os.path.basename(ug_location)
    if ug_model_name.startswith('base_'):
        model_name = ug_model_name[5:]
    elif ug_model_name.startswith('ug_'):
        model_name = ug_model_name[3:]
    else:
        model_name = ug_model_name
    gf_location = os.path.join(args.output_folder,'gapfilled_models','gf_{}'.format(model_name))
    write_sbml_model(cobra_model = gf_model, filename = gf_location)

def create_output_folder(args):
    if os.path.exists(os.path.dirname(args.output_folder)):
        if args.fasta_folder:
            base_model_folder = os.path.join(args.output_folder, 'base_models')
            if not os.path.exists(base_model_folder):
                os.makedirs(base_model_folder, exist_ok=True)
            else:
                print('# WARNING: base models folder allready exists')

        gf_model_folder = os.path.join(args.output_folder, 'gapfilled_models')
        if not os.path.exists(gf_model_folder):
            print("# Creating output_folder for gapfilled models")
            os.makedirs(gf_model_folder, exist_ok=True)
        else:
            print('# WARNING: gapfilled models folder allready exists')
    else:
        sys.exit('# ERROR: your output_folder should have parents')

def main():
    args = parse_arguments()
    lazy_loading_modules()
    create_output_folder(args)

    if args.fasta_folder:
        list_of_genomes = [i for i in os.listdir(args.fasta_folder) if i.endswith('.faa')]
        if len(list_of_genomes) == 0:
            sys.exit('# ERROR: no fasta files found in {}'.format(args.fasta_folder))
        print('# Building and gapfilling {} models'.format(len(list_of_genomes)))
        for genome in list_of_genomes:
            print('# Building and gapfilling: {}'.format(genome))
            build_gapfilled_model_from_fasta(os.path.join(args.fasta_folder, genome), args)
        print('# Done')

    elif args.model_folder:
        list_of_models = [i for i in os.listdir(args.model_folder) if i.endswith('.xml')]
        if len(list_of_models) == 0:
            sys.exit('# ERROR: no xml files found in {}'.format(args.model_folder))
        print('# Gapfilling {} models'.format(len(list_of_models)))
        for model in list_of_models:
            print('# Gapfilling: {}'.format(model))
            gapfill_model_wrapper(os.path.join(args.model_folder, model), args)
        print('# Done')
    else:
        sys.exit('I dont think this message can show, if it does, you did something weird')


if __name__ == '__main__':
    main()
