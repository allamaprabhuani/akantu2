//
//  contact.hh
//  akantu
//
//  Created by Alejandro Marcos Aragón on 6/24/14.
//
//

#ifndef __AKANTU_CONTACT_HH__
#define __AKANTU_CONTACT_HH__


#include "contact_common.hh"
#include "parser.hh"
#include "search.hh"
#include "resolution.hh"

__BEGIN_AKANTU__

template <int Dim, template <int> class Search_policy, class RP>
class Contact : public Search_policy<Dim>, public ContactResolution <Dim, RP::analysis, RP::method> {
  
	typedef SolidMechanicsModel model_type;
	typedef Search_policy<Dim> search_type;
	typedef ContactResolution <Dim, RP::analysis, RP::method> resolution_type;

public:
	Contact(int argc, char *argv[], model_type & m) : search_type(m), resolution_type(m)
	{
    // read parameters from file
		std::pair <Parser::const_section_iterator, Parser::const_section_iterator>
		sub_sect = getStaticParser().getSubSections(_st_contact);

		if (sub_sect.first != sub_sect.second)
			this->parseSection(*sub_sect.first);

    // read parameters the command line
		contact_argparser.parse(argc, argv, cppargparse::_remove_parsed);
    
    // finish initialization of resolution class
    this->initialize();
	}


	//! Provide standard output of contact object
	friend std::ostream& operator << (std::ostream & os, const Contact &cd) {
		const resolution_type& r(cd);
		const search_type& s(cd);

		os << "\nContact object info:" << endl;
		os << "  Search type: " << s << endl;
		os << "  Resolution type: " << r << endl;
		return os;
		}
};


template <ContactImplementationMethod i, class contact_type>
void solveContactStep(contact_type& c)
{ c.solveContactStep<i>(&c); }



__END_AKANTU__


#endif /* __AKANTU_CONTACT_HH__ */
