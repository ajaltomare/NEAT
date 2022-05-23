import click
@click.command()
#@click.argument('filename', type=click.Path(exists=True))
#def touch(filename):
#	"""Print FILENAME if the file exists."""
#	click.echo(click.format_filename(filename))
#@click.option('--R', default = 100, help = 'The desired read length')
#@click.option('--o', default = /home/suvinit/Desktop/command-script-test-run-H1N1-1, help = 'output_prefix')

@click.option('--bam', default = "/home/suvinit/NEAT/testing_command_script/test-output-H1N1.bam", help='output golden BAM file')
def hi(bam):
	click.echo(f'\noutput bam test line:\n {bam}')
if __name__ == "__main__":
	hi()
