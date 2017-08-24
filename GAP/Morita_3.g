LoadPackage("CTblLib");
LoadPackage("edim");

#Tries to identify a block by given charactertable and decomposition matrix
BlockIdentify:=function(dm, nb)
	local i;
	for i in [1..Length(nb)] do
		if Length(dm)=nb[i][3] and Length(dm[1])=nb[i][4] and nb[i][5]=Length(Positions(ElementaryDivisorsMat(TransposedMat(dm)*dm), 1)) then
			if nb[i][6] then return i;fi;
			return -i;
		fi;
	od;
	return 0;
end;

p:=3;
#The Morita equivalence classes of the normal blocks. Format: [GroupID, NrBlock, k(B), l(B), elementary divisors 1, data is unique?]
nblocks:=[
[[9,2], 1, 9, 1, 0, false],
[[18,3], 1, 9, 2, 0, true],
[[18,4], 1, 6, 2, 1, false],
[[36,9], 1, 6, 4, 3, true],
[[36,10], 1, 9, 4, 1, true],
[[72,23], 2, 6, 1, 0, true],
[[72,39], 1, 9, 8, 7, true],
[[72,40], 1, 9, 5, 2, true],
[[144,117], 2, 6, 2, 1, false],
[[72,41], 1, 6, 5, 4, true],
[[144,182], 1, 9, 7, 5, true]
];

#A list of different Morita equivalence classes we find by investigating (modular) character tables.
#An entry consists of [DecompositionMatrix,CartanMatrix Suspected derived equivalence class, LIST], where LIST is a list of groups possibly having a ME block. LIST has entries [GroupID, BlockNr]
mc:=[];
for entry in nblocks do
	ct:=CharacterTable(SmallGroup(entry[1]));
	mct:=ct mod p;
	Q:=DecompositionMatrix(mct, entry[2]);
	C:=TransposedMat(Q)*Q;
	Add(mc, [Q,C, BlockIdentify(Q, nblocks), [ [entry[1], entry[2]] ] ]);
od;

for x in AllCharacterTableNames(IsDuplicateTable, false) do
    ct:=CharacterTable(x);
    if Size(ct) mod p*p=0 and (not IsPSolvableCharacterTable(ct,p)) and ct mod p<>fail then       
	pb:=PrimeBlocks(ct,p);
        mct:=ct mod p;
        for i in [1..Size(pb.defect)] do
		if pb.defect[i]<>2 then break;fi;
		Q:=DecompositionMatrix(mct, i);
		C:=TransposedMat(Q)*Q;
		match:=false;
        	for entry in mc do
			if TransformingPermutations(C, entry[2])<>fail then
				Add(entry[4], [x,i]);
				match:=true;
			fi;
		od;
		if match=false then
			Add(mc, [Q,C, BlockIdentify(Q, nblocks), [[x, i]] ]);
		fi;
        od;
    fi;
od;

#Unsorted output
for entry in mc do
	Print(entry[3]," ",entry[4],"\n");
od;
#Output sorted in suspected DE classes
for i in [1..Length(nblocks)] do
	Print("B",i,"\n");
	counter:=1;
	for entry in mc do
		if entry[3]=i or entry[3]=-i then 
			Print(counter," ", entry[4],"\n");
			counter:=counter+1;
		fi;
	od;
od;
QUIT_GAP();
