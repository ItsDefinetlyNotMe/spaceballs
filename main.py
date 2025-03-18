import click
import os
from SpaceBalls import SpaceBall

@click.command()
@click.argument('mode', type=click.Choice(['editor', 'plot'], case_sensitive=False))
@click.argument('path', type=click.Path(dir_okay=True, readable=True))
@click.option('-p', '--printd', is_flag=True, help="Flag for printing details (used when 'plot' mode is selected)")
@click.option('-e', '--epsilon', type=float, help="time between collisions")
@click.option('-s', '--speed', is_flag=True, help="plot speed")
@click.option('-v', '--velocity', is_flag=True, help="plot velocity")
@click.option('-d', '--distance', is_flag=True, help="plot distance traveled")
@click.option('-r', '--reverse', is_flag=True, help="reverses Plot axis")
@click.option('-my', '--maxy', type=float, help='maximum y displayed')
@click.option('-ly', '--miny', type=float, help='minimum y displayed')
@click.option('-x', '--highlightx', type=float, multiple=True, help='highlight x values')


def main(mode, path, printd, epsilon, speed, velocity, distance, reverse, highlightx, maxy, miny):
    arguments = {'p': True, 's': speed, 'v': velocity, 'd': distance, 'r': reverse, 'hx': highlightx, 'my':maxy, 'ly':miny}
    if epsilon is None:
        epsilon = 0.0001
    VM = SpaceBall(epsilon, arguments)
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