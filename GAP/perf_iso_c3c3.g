PerfIsoTest:=function(ctG, ctH, blG, blH, sinG, sinH, regG, regH, perm, sign)
	local i, j, k, mu;
	#tests separation condition
	#we are testing this condition fist, because it rules out the most cases
	for i in regG do
		for j in sinH do
			mu:=0;		
			for k in [1..Length(blG)] do
				mu:=mu+ValuesOfClassFunction(Irr(ctG)[blG[k]])[i]*sign[k]*ValuesOfClassFunction(Irr(ctH)[blH[perm[k]]])[j];
			od;
			if not mu=0 then return false;fi;
		od;
	od;
	#we even skip the case of considering i in sinG and j in regH as it does not rule out anything further
	#isometry condition
	for i in [1..Length(blG)] do
		for j in [1..Length(blG)] do
			if ScalarProduct(Irr(ctG)[blG[i]], Irr(ctG)[blG[j]])<>sign[i]*sign[j]*ScalarProduct(Irr(ctH)[blH[perm[i]]], Irr(ctH)[blH[perm[j]]]) then Print(perm, " ", sign, " sp", "\n");return false; fi;
		od;
	od;
	#integrality condition
	#since testing on algebraic integer is not quite easily done we will let the given perm/sign pass the test, if it is mu=0, 9, -9 since these are certainly algebraic even if divided by 9 (order of the centralizer). In case it is not we will let the test fail and print mu, testing mu on being an algebraic interger has then to be done by brain. Luckily this is pretty straightforward and only happens in few cases. Just so we are able to check we also print the permutation, signs and conjugacy classes in that case.
	for i in sinG do
		for j in sinH do
			mu:=0;		
			for k in [1..Length(blG)] do
				mu:=mu+ValuesOfClassFunction(Irr(ctG)[blG[k]])[i]*sign[k]*ValuesOfClassFunction(Irr(ctH)[blH[perm[k]]])[j];
			od;
			if mu<>0 and mu<>9 and mu<>-9 then Print(perm, " ", sign, " ", i, " ", j," mu=", mu, "\n");return false;fi;
		od;
	od;
	return true;
end;

#Initialize our situation
G:=SmallGroup(18, 4);	#B3 is the principle block of this group
H:=SmallGroup(144, 117);#b10 is the (unique) non-principle block of this group
ctG:=CharacterTable(G);
ctH:=CharacterTable(H);
pbG:=PrimeBlocks(ctG, 3);
pbH:=PrimeBlocks(ctH, 3);
blG:=Positions(pbG.block, 1);
blH:=Positions(pbH.block, 2);
#Take a look at the character tables for the two blocks
Display(ctG, rec(chars:=blG));
Display(ctH, rec(chars:=blH));

sinG:=[3,4,5,6];	#singular conjugacy classes of G
regG:=[1,2];		#regular conjugacy classes of G
sinH:=[6,8,9,10,11,13,14,15];	#singular conjugacy classes of H
regH:=[1,2,3,4,5,7,12];		#regular conjugacy classes of H

perms:=Arrangements([1,2,3,4,5,6], 6);	#all possible permutations of the characters
signs:=Tuples([1,-1], 6);		#all possible signs for a permutation of characters

#Now check the conditions for every possible permutation/sign combination by using our PerfIsoTest
for perm in perms do
for sign in signs do
	if PerfIsoTest(ctG, ctH, blG, blH, sinG, sinH, regG, regH, perm, sign) then Print(perm," ", sign," ", "\n");fi;
od;
od;
QUIT_GAP();
