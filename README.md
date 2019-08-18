## Usage ##

The package _KhT_ computes the Khovanov and Bar-Natan homology of arbitrary tangles from a splitting into elementary tangles. For pointed 4-ended tangles, it also computes the immersed curve invariants from [[KWZ]](https://github.com/spinachstealer/KhT#References). An example can be found in `examples/2-cable-trefoil.py`. Simply open the terminal in the root folder and run

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

## References ##

[KWZ] Artem Kotelskiy, Liam Watson and Claudius Zibrowius: Immersed curve in Khovanov homology, in preparation. 

## License ##

This program is licensed under the GNU General Public License, version 3. 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
