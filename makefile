all : 
	g++ -flto -Ofast -DNDEBUG -march=native -o DnBk main.cpp Graph.cpp -w
debug:
	g++ -g -O0 -march=native -o DnBk main.cpp Graph.cpp -w
sub:
	g++ -flto -Ofast -DNDEBUG -o MCC SubGraph.cpp
clean:
	rm -rf DnBk
test:
	g++ -flto -Ofast -DNDEBUG -march=native -o DnBk main.cpp Graph.cpp -w
	./DnBk ./data/socfb-Berkeley13.bin 1
	./DnBk ./data/socfb-Texas84.bin 1
	./DnBk ./data/soc-youtube.bin 1