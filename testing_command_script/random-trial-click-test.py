import click
@click.command()
@click.argument('filename', type=click.Path(exists=True))
def touch(filename):
	"""Print FILENAME if the file exists."""
	click.echo(click.format_filename(filename))

@click.option('--R', default = 100, help = 'The desired read length')
@click.option('--o', default = /home/suvinit/Desktop/command-script-test-run-H1N1-1, help = 'output_prefix')
