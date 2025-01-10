#!/usr/bin/env python3
import click


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "--filename", "-f", default="module.txt", help="The name of the file to write to."
)
def write_file(filename):
    with open(filename, "w") as f:
        f.write("Hello, world!\n")


if __name__ == "__main__":
    write_file()
