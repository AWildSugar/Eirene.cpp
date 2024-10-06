vr: 
	g++ -fPIC -fno-exceptions -o ./build/test_vr.out -O0 -g -Wall -Wextra -fsanitize=address -lfmt -std=c++20 src/test_vr.cpp src/eirene.cpp src/util.cpp src/eirene_impl.cpp 
eirene: 
	g++ -fPIC -fno-exceptions -o ./build/eirene_test.out -O0 -g -fsanitize=address -Wall -Wextra -lfmt -std=c++20 src/eireneTest.cpp src/eirene.cpp src/eirene_impl.cpp src/util.cpp -I/opt/homebrew/include/ -L/opt/homebrew/lib

staticlib:
	g++ -fPIC -fno-exceptions -c -o ./build/eirene.o -Wall -Wextra -lfmt -std=c++20 src/eirene.cpp
	ar rc lib/libeirene.a ./build/eirene.o 

sharedlib:
	g++ -fPIC -shared -fno-exceptions -o lib/libeirene.so -Wall -Wextra -lfmt -std=c++20 src/eirene.cpp -I/opt/homebrew/include/


minerva: 
	g++ -fPIC -fno-exceptions -o ./build/minerva_test.out -O0 -g -fsanitize=address -Wall -Wextra -lfmt -std=c++20 src/eirene.cpp src/minerva.cpp src/test_minerva.cpp src/eirene_impl.cpp src/util.cpp -I/opt/homebrew/include/ -L/opt/homebrew/lib
