#!/usr/bin/env python3

import iterm2
import AppKit

# Launch the app
AppKit.NSWorkspace.sharedWorkspace().launchApplication_("iTerm2")

async def main(connection):

    # Get app
    print("Test::main")
    app = await iterm2.async_get_app(connection)

    # Ensure window
    window = app.current_terminal_window
    if app.current_terminal_window is None:
        exit()

    # Create main tab
    main = await window.async_create_tab()
    await main.async_activate()
    await main.async_set_title('CAVSnet DevOps')
    singapura = main.current_session
    await singapura.async_send_text('ssh atks@singapura\n')

    # Split singapura
    shangrila = await singapura.async_split_pane(vertical=True)
    await shangrila.async_send_text('ssh atks@shangri-la\n')

    # Split singapura
    atlantis = await singapura.async_split_pane(vertical=False)
    await atlantis.async_send_text('ssh atks@atlantis\n')

    # Split shangrila
    shambhala = await shangrila.async_split_pane(vertical=False)
    await shambhala.async_send_text('ssh atks@shambhala\n')

    # Split atlantis
    temasek  = await atlantis.async_split_pane(vertical=False)
    await temasek.async_send_text('ssh atks@temasek\n')


iterm2.run_until_complete(main)
