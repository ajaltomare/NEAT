#1>>redo options with click.option
import click
@click.command()
#no capital letters:
@click.option('--file_path', default = "/home/suvinit/NEAT-data/H1N1/H1N1.fa", help = 'Path to reference fasta', type=click.Path(exists=True))

@click.option('--read_length', default = 100, help = 'The desired read length')
#try to test is ^ this is within a certain range (50,300)

@click.option('--output_path', default = "/home/suvinit/NEAT-data/H1N1/H1N1.test-run", help = 'output_prefix')
@click.option('--bam', default = "/home/suvinit/NEAT/testing_command_script/test-output-H1N1.bam", help='output golden BAM file')
################################<TASK BELOW>################################
#try validating that the file directory exists (/home/suvinit)

def hello(file_path, read_length, output_path, bam):
	click.echo(f'Hello {file_path}' + f'\nHello {read_length}' + f'\nHello {output_path}')
	click.echo(f'\noutput bam test line:\n {bam}')
if __name__ == '__main__':
	hello()
