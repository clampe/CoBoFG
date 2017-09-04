#Tests the separation condition for a given perm/signs for the selected characters
Separation:=function(ctG, ctH, blG, blH, sinG, sinH, regG, regH, perm, sign, selection)
	local i, j, k, mu;
	for i in sinG do
		for j in regH do
			mu:=0;		
			for k in [1..Length(selection)] do
				mu:=mu+ValuesOfClassFunction(Irr(ctG)[blG[selection[k]]])[i]*sign[k]*ValuesOfClassFunction(Irr(ctH)[blH[perm[k]]])[j];
			od;
			if not mu=0 then return false;fi;
		od;
	od;
	return true;
end;

#Tests the given polynomial with indeterminates as roots of unity for being an algebraic integer
IntCoeffTest:=function(poly)
	local i, coeff, cshort;
	coeff:=CoefficientsOfUnivariatePolynomial(poly);
	cshort:=[0, 0, 0, 0, 0];
	for i in [1..Length(coeff)] do
		cshort[i mod 5+1]:=cshort[i mod 5+1]+coeff[i];
	od;
	for i in [1..5] do
		if cshort[i] mod 25<>0 then return false;fi;
	od;
	return true;
end;

#Tests if the given perm/sign combination satisfies the integrality condition
Integrality:=function(perm, sign)
	local x, b, c, d, e, f, cb, cc, sd, me, mf, valuesG, valuesH, mu, i, v, w;
	x:=Indeterminate(CF(5), "x");
	b:=2*x+x^3;
	c:=2*x^3+x^4;
	d:=-x-x^4;
	e:=-x^2+x^3;
	f:=-x+x^4;
	cb:=x^2+2*x^4;
	cc:=x+2*x^2;
	sd:=-x^2-x^3;
	me:=x^2-x^3;
	mf:=x-x^4;
	valuesG:=[ [1, 1, 1, b, c, cc, cb, d, d, sd, sd], [1, 1, 1, sd, d, d, sd, cb, b, c, cc], [1, 1, 1, d, sd, sd, d, c, cc, b, cb] ];
	valuesH:=[ [e, me, f, me], [mf, f, e, me], [me, e, mf, f] ];
	for v in valuesG do
		for w in valuesH do
			mu:=0;
			for i in [1..4] do
				mu:=mu+w[i]*sign[i]*v[perm[i]];
			od;
			if not IntCoeffTest(mu) then return false;fi; 
		od;
	od;
	return true;
end;

#Initialize the present situation
H:=SmallGroup(75, 2);
G:=SmallGroup(600, 61);
ctG:=CharacterTable(G);
ctH:=CharacterTable(H);
pbG:=PrimeBlocks(ctG, 5);
pbH:=PrimeBlocks(ctH, 5);
blG:=Positions(pbG.block, 2);
blH:=Positions(pbH.block, 1);

Display(ctH, rec(chars:=blH));
Display(ctG, rec(chars:=blG));

sinH:=[3,4,6,7,8,9,10,11];
regH:=[1,2,5];
sinG:=[19,20,21,22];
regG:=[1,6,7,12,13,18,23,24,25];

#Select some characters with many zero entries to simplyfy calculations
selection:=[8,9,10,11];
Print("Selected characters of G: ");
for i in selection do
	Print(blG[i]," ");
od;
Print("\n");

images:=Combinations([1,2,3,4,5,6,7,8,9,10,11], 4);	#Alle moeglichen Isometrien beschraenkt auf 18,19,20,21 und sortiert ohne Vorzeichen
signs:=Tuples([-1,1],4);				#Alle moeglichen Vorzeichenkombinationen
passed:=[];						#Hier sammel ich alle Permutationen bei denen \mu(g,h)=0 fuer g aus selG und h aus selH(d.h. mit Sicherheit eine alg. Zahl)

#Now we check the conditions. In contrast to the case of C_3 x C_3 checking for algebraic integers is now automated
for image in images do
	for perm in Arrangements(image, Length(image)) do
		for sign in signs do
			if Separation(ctG, ctH, blG, blH, sinG, sinH, regG, regH, perm, sign, selection) and Integrality(perm, sign) then Add(passed, [perm, sign]);fi;
		od;
	od;
od;
Print("List of perm/signs that satisfy the separation as well as the Integrality condition:\n");
Print(passed,"\n");
Print("If it is empty, then there are none.\n");
QUIT_GAP();
