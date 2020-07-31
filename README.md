# JULIASCII

A multithreaded ASCII/ANSI Julia fractal viewer with PPM export.

* Somewhat glitchy multithreading when closely inspecting exported PPMs.
* Certain character "pallettes" may not be friendly to UTF-8 terminals.
* Single precision by default, a #define provided for easy configuration.
* Supports "script" execution, but you could as well pipe it in.
* Controls are unlike anything, but an ASCII help chart is included.
* Uses ANSI control sequences for redrawing, NTVDM needed on Windows <10.
* No unbuffered input, but you can "charge" the Return/Enter/<-' key.

This program was originally written in 2018Q2 for a university semester work in Computer Architecture course on some custom Xilinx development board (ARM32), however all functionality specific to this device has been removed, as I no longer have access to it. Being written in C that I started to "git gud" in, it sucked me in so much that I gladly sacrificed the paralell Java course which contained another semester work.

Clearly not as good as XaoS, but at least can be run in a terminal without dependence on either aalib, ncurses, or libcaca.

Generated Doxygen documentation not included, just read the source file, it's not even 700 lines long. Even though, it's still too long to paste in a blog post.

As for a loicence for dat, any version of GPL or CC BY-SA will do. They may be incompatible in pedantic details, but in spirit they are, so take that as double-licencing. Binaries are not provided, so they can't come with any licence.
