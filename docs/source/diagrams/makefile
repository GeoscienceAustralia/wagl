

pq.png:	pq.dot
	dot -Tpng -l DatabaseShape.ps pq.dot -o pq.png
pq.ps:	pq.dot
	dot -Tps -l DatabaseShape.ps pq.dot -o pq.ps

pq.pdf:	pq.ps
	estopdf pq.ps

all:	pq.png pq.ps pq.pdf
