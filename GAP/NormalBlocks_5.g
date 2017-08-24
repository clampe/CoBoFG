Read("ListNormalBlocks.g");

#Print the blocks we obtained with their inertial quotients, group, schur covers, k(B) and l(B)
p:=5;
Print("List for p=",p,": \n");
Print("Blocks with defect group C_p \n");
ListCyclicBlocks(p, 1);
Print("\nBlocks with defect group C_p^2 \n");
ListCyclicBlocks(p, 2);
Print("\nBlocks with defect group C_p x C_p\n");
ListNormalBlocks(p);
Print("\n");


#The next calculation is needed for cases 3 and 4 in the proof of the normal classification
Print("Calculations for cases 3 and 4 in the normal classification proof: \n\n");
#We need to show that all automorphisms of Z are restrictions of automorphisms of the schur cover S
#Case 3
G:=SmallGroup(1600, 5606);
C:=Centre(G);
gen:=Pcgs(C);
imC:=[];
imG:=[];
Print(StructureDescription(C), " is cyclic  with generator ", gen[1], " \n");
for aut in AutomorphismGroup(C) do
	AddSet(imC, Image(aut, gen[1]));
od;
for aut in AutomorphismGroup(G) do
	if Length(imG)=Length(imC) then break;fi;
	AddSet(imG, Image(aut, gen[1]));
od;
Print("Images of the generator under automorphisms of Z: \n");
for entry in imC do
	Print(entry, ", ");
od;
Print("\n");
Print("Images of the generator under automorphisms of S: \n");
for entry in imG do
	Print(entry, ", ");
od;
Print("\nSince they are equal, we see that all automorphisms of Z are restrictions of automorphisms of S\n\n");

#Case 4
G:=SmallGroup(1600, 5725);
C:=Centre(G);
gen:=Pcgs(C);
imC:=[];
imG:=[];
Print("Generators of ", StructureDescription(C), " :", gen[1]," ", gen[2], " \n");
for aut in AutomorphismGroup(C) do
	AddSet(imC, [Image(aut, gen[1]), Image(aut, gen[2])]);
od;
for aut in AutomorphismGroup(G) do
	if Length(imG)=Length(imC) then break;fi;
	AddSet(imG, [Image(aut, gen[1]), Image(aut, gen[2])]);
od;
Print("Images of the generators under automorphisms of Z: \n");
for entry in imC do
	Print(entry, ", ");
od;
Print("\n");
Print("Images of the generators under automorphisms of S: \n");
for entry in imG do
	Print(entry, ", ");
od;
Print("\nSince they are equal, we see that all automorphisms of Z are restrictions of automorphisms of S\n\n");
#QUIT_GAP();
