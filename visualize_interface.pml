load 6m0j_fixed.pdb
bg_color white
hide all
show cartoon
color gray80
select interface_A, chain A and resi 19+24+27+28+30+31+34+35+37+38+41+42+45+79+82+83+324+325+326+330+353+354+355+357+393
show sticks, interface_A
color cyan, interface_A
color deepteal, interface_A and name C*
select interface_E, chain E and resi 403+417+445+446+449+453+455+456+473+475+476+477+484+486+487+489+493+496+498+500+501+502+503+505
show sticks, interface_E
color orange, interface_E
color brightorange, interface_E and name C*
deselect
zoom interface_A or interface_E, 5
set transparency, 0.4
