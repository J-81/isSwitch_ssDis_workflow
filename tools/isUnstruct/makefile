#include "3_1.hpp"
#include "pdb_data.hpp"
#include "vector.hpp"
#include "str_data.hpp"
#include "Phy-Psy.hpp"
#include "my-cgi.hpp"
#include "ReadArray.hpp"
#include "rnd.hpp"
#include "stat.hpp"
#include "min_max.hpp"
#include "put_exp.hpp"
#include "buffer.hpp"

CDIR= /home/mlobanov/cprog/COmON/
CFLAGS= -c -O3 -I/home/mlobanov/cprog/COmON/
OBJ_GL= $(CDIR)pdb_data.o $(CDIR)str_data.o $(CDIR)vector.o $(CDIR)3_1.o  $(CDIR)Phy-Psy.o $(CDIR)my-cgi.o $(CDIR)ReadArray.o $(CDIR)rnd.o $(CDIR)stat.o $(CDIR)min_max.o $(CDIR)put_exp.o $(CDIR)buffer.o
GL_HPP= $(CDIR)pdb_data.hpp $(CDIR)str_data.hpp $(CDIR)vector.hpp $(CDIR)3_1.hpp $(CDIR)Phy-Psy.hpp $(CDIR)my-cgi.hpp $(CDIR)ReadArray.hpp $(CDIR)rnd.hpp $(CDIR)stat.hpp $(CDIR)min_max.hpp $(CDIR)put_exp.hpp $(CDIR)buffer.hpp
IsUnstruct : pred_del.o is_unstruct.o
	g++ -lm -O3 -o IsUnstruct pred_del.o is_unstruct.o
pred_del.o : pred_del.cpp pred_del.hpp
	g++ $(CFLAGS) pred_del.cpp
is_unstruct.o : is_unstruct.cpp pred_del.hpp
	g++ $(CFLAGS) is_unstruct.cpp

$(CDIR)3_1.o : $(CDIR)3_1.cpp $(CDIR)3_1.hpp
	g++ $(CFLAGS) -o $(CDIR)3_1.o $(CDIR)3_1.cpp
$(CDIR)pdb_data.o : $(CDIR)pdb_data.cpp $(CDIR)pdb_data.hpp
	g++ $(CFLAGS) -o $(CDIR)pdb_data.o $(CDIR)pdb_data.cpp
$(CDIR)vector.o : $(CDIR)vector.cpp $(CDIR)vector.hpp
	g++ $(CFLAGS) -o $(CDIR)vector.o $(CDIR)vector.cpp
$(CDIR)str_data.o : $(CDIR)str_data.cpp $(CDIR)str_data.hpp
	g++ $(CFLAGS) -o $(CDIR)str_data.o $(CDIR)str_data.cpp
$(CDIR)Phy-Psy.o : $(CDIR)Phy-Psy.cpp $(CDIR)Phy-Psy.hpp
	g++ $(CFLAGS) -o $(CDIR)Phy-Psy.o $(CDIR)Phy-Psy.cpp
$(CDIR)my-cgi.o : $(CDIR)my-cgi.cpp $(CDIR)my-cgi.hpp
	g++ $(CFLAGS) -o $(CDIR)my-cgi.o $(CDIR)my-cgi.cpp
$(CDIR)ReadArray.o : $(CDIR)ReadArray.cpp $(CDIR)ReadArray.hpp
	g++ $(CFLAGS) -o $(CDIR)ReadArray.o $(CDIR)ReadArray.cpp
$(CDIR)rnd.o : $(CDIR)rnd.cpp $(CDIR)rnd.hpp
	g++ $(CFLAGS) -o $(CDIR)rnd.o $(CDIR)rnd.cpp
$(CDIR)stat.o : $(CDIR)stat.cpp $(CDIR)stat.hpp
	g++ $(CFLAGS) -o $(CDIR)stat.o $(CDIR)stat.cpp
$(CDIR)min_max.o : $(CDIR)min_max.cpp $(CDIR)min_max.hpp
	g++ $(CFLAGS) -o $(CDIR)min_max.o $(CDIR)min_max.cpp
$(CDIR)put_exp.o : $(CDIR)put_exp.cpp $(CDIR)put_exp.hpp
	g++ $(CFLAGS) -o $(CDIR)put_exp.o $(CDIR)put_exp.cpp
$(CDIR)buffer.o : $(CDIR)buffer.cpp $(CDIR)buffer.hpp
	g++ $(CFLAGS) -o $(CDIR)buffer.o $(CDIR)buffer.cpp
