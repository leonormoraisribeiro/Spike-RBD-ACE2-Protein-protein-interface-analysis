load 6m0j_fixed.pdb
bg_color white
hide all
show cartoon
color gray90
set transparency, 0.6
select core_interface, (chain A and resi 19+24+27+28+30+31+324+325+326+330+34+35+353+354+355+357+37+38+386+393+41+42+45+79+82+83) or (chain E and resi 403+417+445+446+449+453+455+456+473+475+476+477+484+485+486+487+489+490+493+496+498+500+501+502+503+505)
show sticks, core_interface
color magenta, core_interface
color hotpink, core_interface and name C*
select periphery_dist, (chain A and resi 101+20+21+23+25+26+29+32+327+329+33+331+351+352+356+36+383+387+388+389+39+390+40+43+44+46+48+49+75+76+78+80+81+84+97) or (chain E and resi 405+406+408+409+418+421+439+443+444+447+448+454+457+474+478+488+491+492+494+495+497+499+504+506+507)
show lines, periphery_dist
color lightblue, periphery_dist
deselect
zoom core_interface, 5
set ray_shadows, 0
set stick_radius, 0.25
echo [LEGEND] Magenta: Core Interface (Energy + Dist)
echo [LEGEND] LightBlue: Periphery (Dist Only)
