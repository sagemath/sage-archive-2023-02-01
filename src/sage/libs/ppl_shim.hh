#ifndef PPL_SHIM__H
#define PPL_SHIM__H



#include <ppl.hh>

using namespace Parma_Polyhedra_Library;


// access Generator's static methods
Generator* new_line(const Linear_Expression& e);
Generator* new_ray(const Linear_Expression& e);
Generator* new_point(const Linear_Expression& e, Coefficient d);
Generator* new_closure_point(const Linear_Expression& e, Coefficient d);

// Poly_Gen_Relation/Poly_Con_Relation have no default constructor
Poly_Gen_Relation* new_relation_with(const Polyhedron &p, const Generator &g);
Poly_Con_Relation* new_relation_with(const Polyhedron &p, const Constraint &c);


// Iterator for Generator_System
typedef Generator_System::const_iterator* gs_iterator_ptr;
gs_iterator_ptr init_gs_iterator(const Generator_System &gs);
Generator next_gs_iterator(gs_iterator_ptr);
bool is_end_gs_iterator(const Generator_System &gs, gs_iterator_ptr gsi_ptr);
void delete_gs_iterator(gs_iterator_ptr);


// Iterator for Constraint_System
typedef Constraint_System::const_iterator* cs_iterator_ptr;
cs_iterator_ptr init_cs_iterator(const Constraint_System &cs);
Constraint next_cs_iterator(cs_iterator_ptr);
bool is_end_cs_iterator(const Constraint_System &cs, cs_iterator_ptr csi_ptr);
void delete_cs_iterator(cs_iterator_ptr);



#endif
