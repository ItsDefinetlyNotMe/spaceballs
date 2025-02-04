import click
import os
from SpaceBalls import SpaceBall

@click.command()
@click.argument('mode', type=click.Choice(['editor', 'plot'], case_sensitive=False))
@click.argument('path', type=click.Path(dir_okay=True, readable=True))
@click.option('-p', '--printd', is_flag=True, help="Flag for printing details (used when 'plot' mode is selected)")
@click.option('-e', '--epsilon', type=float, help="time between collisions")
def main(mode, path, printd, epsilon):
    if epsilon is None:
        epsilon = 0.0001
    VM = SpaceBall(epsilon)
    if mode == 'editor':
        if not os.path.exists(path):
            click.echo("Path not found opening empty File.")
            VM.add_live_scene()
        else:
            VM.add_live_scene(path)
    elif mode == 'plot':
        if not os.path.exists(path):
            click.echo("Path does not exist.")
            exit()
        else:
            VM.add_scene(path)
            VM.play_sequences()
            VM.plot_plot()

            if printd:
                VM.print_details()
            
            VM.show()



if __name__ == '__main__':
    main()