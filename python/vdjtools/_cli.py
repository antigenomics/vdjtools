"""Set signature kernel limits before importing the command implementation."""
def main():
    from .cores import configure_cli_runtime
    import sys
    from typer.main import get_command
    from .cli import app

    command = get_command(app).commands.get(sys.argv[1]) if len(sys.argv) > 1 else None
    if command is not None and any(p.name == 'threads' for p in command.params):
        configure_cli_runtime()
    app()
