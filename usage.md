An example computation can be found in `examples/2-cable-trefoil.py`. Simply open the terminal in the root folder and run

	python3 KhT 2-cable-trefoil

This creates the two files

* `examples/2-cable-trefoil_BNr_field=7.pdf` and 
* `examples/2-cable-trefoil_Khr_field=7.pdf`

which contain the reduced Bar-Natan homology and the reduced Khovanov homology of the tangle, respectively.

To compute a new example, make a copy of the template file

	cp examples/template.py examples/new_tangle.py

and edit the new file appropriately. The computation for the new tangle can be run via

	python3 KhT new_tangle

If you have comments or questions about this package, [please get in touch](https://cbz20.raspberryip.com/).
