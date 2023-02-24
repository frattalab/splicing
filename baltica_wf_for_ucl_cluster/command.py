#!/usr/bin/env python
"""
Command line interface for Baltica

"""
import os
from pathlib import Path
import sys
import yaml
import logging
import subprocess
import tempfile
import snakemake
import click

from baltica.version import _program, __version__

baltica_path = Path(__file__)


# from https://stackoverflow.com/a/56944256/1694714
class CustomFormatter(logging.Formatter):
    """Logging Formatter to add colors and count warning / errors"""

    grey = "\x1b[38;21m"
    green = "\x1b[32;21m"
    yellow = "\x1b[33;21m"
    red = "\x1b[31;21m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "baltica:\t| %(message)s"

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: green + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


def avaiable_workflows(baltica_path):
    smk_files = baltica_path.parent.glob('*.smk')
    return [x.with_suffix('').name for x in smk_files]


avaiable_workflows_ = avaiable_workflows(baltica_path)

logger = logging.getLogger(__file__)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(CustomFormatter())
logger.addHandler(ch)

# https://click.palletsprojects.com/en/8.0.x/advanced/#forwarding-unknown-options
# unknown options are passed to snakemake


@click.command(context_settings=dict(ignore_unknown_options=True))
@click.version_option(__version__, prog_name='baltica')
@click.argument("workflow", type=click.Choice(avaiable_workflows_, case_sensitive=False))
@click.argument("config_file",  type=click.Path(exists=True))
@click.option('-v', '--verbose', is_flag=True, help='Enables verbose mode.')
@click.argument('snakemake_args', nargs=-1, type=click.UNPROCESSED)
def cli(workflow, config_file, verbose, snakemake_args):
    """
    Baltica implements workflows for differential junction
    usage methods, and method integration and analysis. Visit
    https://github.com/dieterich-lab/Baltica for more information.


    Runs baltica <WORKFLOW> with <CONFIG_FILE> and <SNAKEMAKE_ARGS>
    """
    # TODO add link to baltica docs with important snakemake parameters
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.ERROR)
    # Error if older version of snakemake is installed
    min_snakemake_version = "6"
    try:
        snakemake.utils.min_version(min_snakemake_version)
    except snakemake.exceptions.WorkflowError:
        logger.error(
            f'{_program} requires Snakemake version >= {min_snakemake_version}:',
            exc_info=True)
        sys.exit(1)
    # check if workflow file is readable
    snakefile = (baltica_path.parent / workflow).with_suffix(".smk")

    with open(config_file) as fin:
        config = yaml.safe_load(fin)

    config['config_path'] = str(Path(config_file).resolve())
    config['baltica_path'] = str(baltica_path.resolve())

    logger.debug(f"Config file is {config['config_path']}")
    with open(config_file, 'w') as fou:
        yaml.dump(config, fou)

    target = config["path"]

    try:
        os.makedirs(Path(target) / 'logs/')
    except FileExistsError:
        pass

    snakemake_args = list(snakemake_args)

    if verbose:
        snakemake_args.extend(['--printshellcmds', '--verbose', '--reason'])
    else:
        logger.warning(
            "Starting Baltica run in a quiet mode. Use --verbose "
            "to change this behavior.")
    if not any([x in snakemake_args for x in ['--cores', '-c', '--job', '-j', '--profile']]):
        logger.warning(
            "Snakemake invoked with a single-core, use --cores N or "
            "--jobs N, where N is the number of available cores to "
            "change this parameter.")
        snakemake_args.append('-j1')

    # Singularity support
    # here we set bindings three directories needed by singularity
    # the target path, where the output is written to
    # the sample path, which contains the input data
    # the baltica directory, which contains the analysis scripts
    if '--use-singularity' in snakemake_args and "--singularity-args" not in snakemake_args:
        relative_path = Path(baltica_path).parent.resolve()
        bound_path = [
            config["path"],
            str(config["sample_path"]),
            str(Path(config['ref']).parent),
            str(Path(config['ref_fa']).parent),
            str(Path(config['config_path']).parent),
            str(relative_path),
            tempfile.gettempdir()]

        if 'orthogonal_result' in config:
            bound_path.append(str(Path(config['orthogonal_result']).parent))

        if 'bind_paths' in config:
            if not isinstance(config['bind_paths'], list):
                logger.error("bind_paths must be a list")
                pass
            bound_path.extend(config['bind_paths'])

        bound_path = set(bound_path)

        # bind several paths that contain input data
        snakemake_args.extend(
            ['--singularity-args', '-B ' + ','.join(bound_path.difference('.'))])

    try:
        _ = subprocess.run(['singularity', '--version'],
                           stdout=subprocess.DEVNULL)
    except FileNotFoundError:
        if '--use-singularity' in snakemake_args:
            logger.critical(
                "Baltica requires Singularity, which was not found", exc_info=True)
            sys.exit(1)

    if '--use-singularity' in snakemake_args and '--singularity-prefix' not in snakemake_args:
        # set $HOME/.baltica/singularity/ as download directory for the containers
        snakemake_args.extend(
            ['--singularity-prefix',
                str(Path.home() / '.baltica/singularity/')]
        )

    if workflow == 'all':
        # append final rule name for end-to-end execution
        logger.info("Running baltica in the end-to-end mode.")
        snakemake_args.append('final')

    logger.info(
        f"""Starting baltica (v{__version__}) analysis with:
    WORKFLOW: {workflow} (from {snakefile})
    CONFIGURATION: {config_file}
    TARGET DIRECTORY: {target}
    SNAKEMAKE ARGUMENTS: {' '.join(snakemake_args)}
    """)
    cmd = [
        'snakemake',
        '--configfile', config_file,
        '--snakefile', str(snakefile),
        *snakemake_args]
    logger.debug('Start of snakemake logger:')
    result = subprocess.run(cmd)
    return result.check_returncode


if __name__ == '__main__':
    cli()
