#!/usr/bin/env python3

__author__ = 'M.D. Boer'
__version__ = '0.1'
__date__ = '15 Jan, 2024'

import argparse
import os
import sys
import logging
import logging.config

def parse_arguments():
    class PathAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            path = os.path.expanduser(values.rstrip('/'))

            if not path.startswith('/') and not path.startswith('.'):
                path = f'./{path}'

            setattr(namespace, self.dest, path)

    parser = argparse.ArgumentParser(
            prog='fasta2model',
            description=('Command line script to build models from a folder with fasta files (-f DIR) or '
            'a folder with ungapfilled base models (-m DIR)'),
            usage=('python fasta2model_CLI.py (-f [DIR] | -d [DIR]) -o output_folder '),
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
                'Folder containing the fasta files you want to build models for'
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

    required.add_argument(
            '-sf',
            '--suffix_faa',
            default='.faa',
            type=str,
            action=PathAction,
            help=('Suffix of the protein fasta')
    )
    
    optional.add_argument(
            '-e',
            '--medium',
            default=None,    
            type=str,
            action=PathAction,
            help=('[File] path to medium file')
    )
    
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

def build_base_model(args):
    from modelseedpy import MSBuilder, MSGenome
    from modelseedpy.core import msmedia
    from modelseedpy.core.rast_client import RastClient
    # Set the path to your genome
    logging.info('# Building MSGenome object')
    patric_genome = MSGenome.from_fasta(args.path_to_fasta, split = '|')
    rast = RastClient()
    rast.annotate_genome(patric_genome)
    logging.info('# Building base model')
    base_model = MSBuilder.build_metabolic_model(model_id = args.model_name,
                                                 genome   = patric_genome,
                                                 index    = '0',
                                                 classic_biomass = True,
                                                 gapfill_model   = False,
                                                 gapfill_media   = None,
                                                 annotate_with_rast = True,
                                                 allow_all_non_grp_reactions = True
                                            )

    logging.info(f'# Base model created with {len(base_model.reactions.list_attr("id"))} reactions')
    return base_model

def build_gapfilled_model_from_fasta(args):
    from cobra.io import write_sbml_model
    args.path_to_base_model = os.path.join(args.output_folder,f'base_models','base_{args.model_name}.xml')
    if not os.path.isfile(args.path_to_base_model):
        base_model = build_base_model(args)
        write_sbml_model(cobra_model = base_model, filename = args.path_to_base_model)
    else:
        logging.warning(f'# Base model {args.model_name} allready exists, skipping')
    gapfill_model_wrapper(args)

def gapfill_model_wrapper(args):
    from cobra.io import write_sbml_model
    from dnngior.gapfill_class import Gapfill
    from dnngior.reaction_class import Reaction

    args.path_to_gf_model = os.path.join(args.output_folder,'gapfilled_models',f'gf_{args.model_name}.xml')
    if not os.path.isfile(args.path_to_gf_model):
        gf_model = Gapfill(args.path_to_base_model, medium_file = args.medium)
        write_sbml_model(cobra_model = gf_model.gapfilledModel, filename =  args.path_to_gf_model)
        n_dr = len(gf_model.draft_reaction_ids)
        n_gr = len(gf_model.added_reactions)
        args.gf_data_file.write(f'{args.model_name}\t{n_dr}\t{n_gr}\t{n_dr+n_gr}\n')
    else:
        logging.warning(f'# Gapfilled model {args.model_name} allready exists, skipping')
        # gf_model = Gapfill(draftModel=args.path_to_gf_model, gapfill=False)
        # gf_model.gapfilledModel = gf_model.draftModel


def create_output_folder(args):
    if os.path.exists(os.path.dirname(args.output_folder)):
        if args.fasta_folder:
            base_model_folder = os.path.join(args.output_folder, 'base_models')
            if not os.path.exists(base_model_folder):
                os.makedirs(base_model_folder, exist_ok=True)
            else:
                print(f'# WARNING: base models folder allready exists ({base_model_folder})')

        gf_model_folder = os.path.join(args.output_folder, 'gapfilled_models')
        if not os.path.exists(gf_model_folder):
            print('# Creating output_folder for gapfilled models')
            os.makedirs(gf_model_folder, exist_ok=True)
        else:
            print(f'# WARNING: gapfilled models folder allready exists ({gf_model_folder})')
        # extra_folder = os.path.join(args.output_folder, 'tmp')
        # if not os.path.exists(extra_folder):
        #     print('# Creating gapfill info folder')
        #     os.makedirs(extra_folder, exist_ok=True)
        # else:
        #     print('# WARNING: tmp folder already exists')
    else:
        sys.exit('# ERROR: your output_folder should have parents')

def main():
    args = parse_arguments()
    create_output_folder(args)
    logging.basicConfig(filename=os.path.join(args.output_folder, 'gapfill.log'), filemode='a', format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

    #I am making my life way too complicated
    gf_data_location = os.path.join(args.output_folder, 'gf_data.tsv')
    if os.path.isfile(gf_data_location):
        file_count = 0
        gf_data_location = os.path.join(args.output_folder, f'gf_data_{file_count}.tsv')
        while os.path.isfile(gf_data_location):
            file_count += 1
            gf_data_location = gf_data_location.replace(str(file_count-1),str(file_count))

    args.gf_data_file = open(gf_data_location, 'w')
    args.gf_data_file.write('model_id\tn_reactions_base\tn_gf_reactions\tn_total_reactions\n')


    if args.fasta_folder:
        list_of_genomes = [i for i in os.listdir(args.fasta_folder) if i.endswith('.faa') or i.endswith('.faa.gz')]
        len_list_of_genomes = len(list_of_genomes)
        if len_list_of_genomes == 0:
            sys.exit(f'# ERROR: no fasta files found in {args.fasta_folder}')
        logging.info(f'# Building and gapfilling {len_list_of_genomes} models')
        for i, genome in enumerate(list_of_genomes):
            args.path_to_fasta = os.path.join(args.fasta_folder, genome)
            args.model_name = '.'.join(os.path.basename(args.path_to_fasta).split('.')[:-1])
            logging.info(f'# Building and gapfilling: {genome} ({i+1}/{len_list_of_genomes})')
            build_gapfilled_model_from_fasta(args)
        logging.info('#Done')

    elif args.model_folder:
        list_of_models = [i for i in os.listdir(args.model_folder) if i.endswith('.xml')]
        len_list_of_models = len(list_of_models)
        if len_list_of_models == 0:
            sys.exit(f'# ERROR: no xml files found in {args.model_folder}')
        logging.info(f'# Gapfilling {len_list_of_models} models')
        for i, model in enumerate(list_of_models):
            args.path_to_base_model = os.path.join(args.model_folder, model)
            args.model_name = '.'.join(os.path.basename(args.path_to_base_model).split('.')[:-1])
            if args.model_name.startswith('base_'):
                args.model_name = args.model_name[5:]
            logging.info(f'# Gapfilling: {model} ({i+1}/{len_list_of_models})')
            gapfill_model_wrapper(args)
        print('# Done')
    else:
        sys.exit('I dont think this message can show, if it does, you (or more likely me) did something weird')
    #args.gf_data_file.close()

if __name__ == '__main__':
    main()
