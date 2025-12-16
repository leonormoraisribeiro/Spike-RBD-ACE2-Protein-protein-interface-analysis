load 6m0j_fixed.pdb
bg_color white
hide all
show cartoon
color gray90
set transparency, 0.5
select dist_A, chain A and resi 101+19+20+21+23+24+25+26+27+28+29+30+31+32+324+325+326+327+329+33+330+331+34+35+351+352+353+354+355+356+357+36+37+38+383+386+387+388+389+39+390+393+40+41+42+43+44+45+46+48+49+75+76+78+79+80+81+82+83+84+97
show sticks, dist_A
color slate, dist_A
select dist_E, chain E and resi 403+405+406+408+409+417+418+421+439+443+444+445+446+447+448+449+453+454+455+456+457+473+474+475+476+477+478+484+485+486+487+488+489+490+491+492+493+494+495+496+497+498+499+500+501+502+503+504+505+506+507
show sticks, dist_E
color raspberry, dist_E
deselect
zoom dist_A or dist_E, 5
echo Geometric Interface (Dist < 8.0 A)
