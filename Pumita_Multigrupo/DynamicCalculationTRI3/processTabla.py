inF= open('TCOMUN_AMPLIADA_original5grupos.txt','r');
outF = open('TCOMUN_AMPLIADA.txt','w')
count = 0;
for line in inF:
	if (len(line)<5):
		count = -2;
	count = count +1;
	if ( not ( (count % 8) in (4,5)  )):
		words=line.split();
		if len(words)>3:
			outF.write(" ".join(words[:3])+'\n');
		else:
			outF.write(line);