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

BlckId:=function(Q)
	local B;
	B:="BX";
	if Size(Q)=13 and Size(Q[1])=10 then B:=28;fi;
	if Size(Q)=20 and Size(Q[1])=14 then B:=36;fi;
	return B;
end;

p:=5;
#The Morita equivalence classes of the normal blocks. Format: [GroupID, NrBlock, k(B), l(B), elementary divisors 1, data is unique?]
nblocks:=[
[[25,2], 1, 25, 1, 0, false],
[[50,3], 1, 20, 2, 0, true],
[[50,4], 1, 14, 2, 1, true],
[[75,2], 1, 11, 3, 2, false],
[[100,9], 1, 25, 4, 0, true],
[[100,10], 1, 13, 4, 2, true],
[[100,11], 1, 10, 4, 3, false],
[[100,12], 1, 10, 4, 3, false],
[[100,13], 1, 16, 4, 1, true],
[[200,25], 2, 13, 1, 0, true],
[[150,5], 1, 13, 3, 1, true],
[[150,6], 1, 10, 6, 5, false],
[[200,40], 1, 11, 8, 7, true],
[[200,41], 1, 20, 8, 3, true],
[[400,118], 2, 14, 2, 0, true],
[[200,42], 1, 14, 8, 5, true],
[[400,125], 2, 8, 2, 1, true ],
[[200,43], 1, 14, 5, 2, true],
[[400,131], 2, 11, 2, 1, true],
[[200,44], 1, 8, 5, 4, true],
[[300,23], 1, 8, 6, 5, true],
[[300,24], 1, 14, 12, 11, true],
[[300,25], 1, 14, 6, 3, true],
[[600,61], 2, 11, 3, 2, false],
[[400,205], 1, 25, 16, 9, true],
[[800,957], 2, 13, 4, 1, true],
[[1600,5606], 3, 10, 1, 0, true],
[[400,206], 1, 13, 10, 8, true],
[[400,207], 1, 16, 10, 6, true],
[[800,968], 2, 10, 4, 2, true],
[[600,148], 1, 13, 12, 11, true],
[[600,149], 1, 25, 24, 23, true],
[[600,150], 1, 8, 7, 6, true],
[[600,151], 1, 16, 12, 9, true],
[[1200,491], 2, 10, 6, 5, false],
[[800,1191], 1, 20, 14, 9, true],
[[1600,9794], 2, 11, 5, 3, true],
[[1200,946], 1, 20, 18, 16, true],
[[1200,947], 1, 16, 14, 12, true],
["PG2519", 1, 20, 16, 12, true]
];

#A list of different Morita equivalence classes we find by investigating (modular) character tables.
#An entry consists of [DecompositionMatrix, Suspected derived equivalence class, LIST], where LIST is a list of groups possibly having a ME block. LIST has entries [GroupID, BlockNr]
mc:=[];
for entry in nblocks do
	if entry[1]="PG2519" then break;fi;
	ct:=CharacterTable(SmallGroup(entry[1]));
	mct:=ct mod p;
	Q:=DecompositionMatrix(mct, entry[2]);
	C:=TransposedMat(Q)*Q;
	Add(mc, [Q,C, BlockIdentify(Q, nblocks), [ [entry[1], entry[2]] ] ]);
od;
ct:=CharacterTable(PrimitiveGroup(25,19));
mct:=ct mod p;
Q:=DecompositionMatrix(mct, 1);
C:=TransposedMat(Q)*Q;
Add(mc, [Q,C, BlockIdentify(Q, nblocks), [ ["PG2519", 1] ] ]);

#Consider all symmetric and alternating groups up to n=18. For bigger groups modular character tables are not available.
for n in [2*p..18] do
	ctS:=CharacterTable(ReplacedString("SB", "B", String(n)));
	ctA:=CharacterTable(ReplacedString("AB", "B", String(n)));
	pbS:=PrimeBlocks(ctS, p);
	pbA:=PrimeBlocks(ctA, p);	
	mctS:=ctS mod p;
	mctA:=ctA mod p;
	#Print("S",n," ",pbS.defect, "\n");
	#Print("A",n," ",pbA.defect, "\n");
	for i in [1..Length(pbS.defect)] do
		if pbS.defect[i]<>2 then continue;fi;
		Q:=DecompositionMatrix(mctS, i);
		C:=TransposedMat(Q)*Q;
		match:=false;
		for entry in mc do
			if TransformingPermutations(C, entry[2])<>fail then
				Add(entry[4], [ReplacedString("S_B", "B", String(n)),i]);
				match:=true;
			fi;
		od;
		if match=false then
			Add(mc, [Q, C, BlckId(Q), [[ReplacedString("S_B", "B", String(n)), i]] ]);
		fi;
	od;
	for i in [1..Length(pbA.defect)] do
		if pbA.defect[i]<>2 then continue;fi;
		Q:=DecompositionMatrix(mctA, i);
		C:=TransposedMat(Q)*Q;
		match:=false;
		for entry in mc do
			if TransformingPermutations(C, entry[2])<>fail then
				Add(entry[4], [ReplacedString("A_B", "B", String(n)),i]);
				match:=true;
			fi;
		od;
		if match=false then
			Add(mc, [Q, C, BlckId(Q), [[ReplacedString("A_B", "B", String(n)), i]] ]);
		fi;
	od;
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
