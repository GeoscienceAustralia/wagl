fc.png:	fc.ps
	convert fc.ps fc.png

fc.ps:	fc.dot
	dot -Tps -l DatabaseShape.ps fc.dot -o fc.ps

fc.pdf:	fc.ps
	epstopdf fc.ps

nbar.png:	nbar.ps
	convert nbar.ps nbar.png

nbar.ps:	nbar.dot
	dot -Tps -l DatabaseShape.ps nbar.dot -o nbar.ps

nbar.pdf:	nbar.ps
	epstopdf nbar.ps

pq.png:	pq.ps
	convert pq.ps pq.png

pq.ps:	pq.dot
	dot -Tps -l DatabaseShape.ps pq.dot -o pq.ps

pq.pdf:	pq.ps
	epstopdf pq.ps

all:	pq.png pq.ps pq.pdf fc.png fc.ps fc.pdf nbar.png nbar.ps nbar.pdf
