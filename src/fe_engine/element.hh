#ifndef __AKANTU_ELEMENT_HH__
#define __AKANTU_ELEMENT_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* Element                                                                    */
/* -------------------------------------------------------------------------- */
class Element;
extern const Element ElementNull;

class Element {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  explicit Element(ElementType type = _not_defined, UInt element = 0,
          GhostType ghost_type = _not_ghost, ElementKind kind = _ek_regular) :
    type(type), element(element),
    ghost_type(ghost_type), kind(kind) {};

  Element(const Element & element) {
    this->type    = element.type;
    this->element = element.element;
    this->ghost_type = element.ghost_type;
    this->kind = element.kind;
  }

  virtual ~Element() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:

  inline bool operator==(const Element & elem) const {
    return ((element == elem.element)
            && (type == elem.type)
            && (ghost_type == elem.ghost_type)
            && (kind == elem.kind));
  }

  inline bool operator!=(const Element & elem) const {
    return ((element != elem.element)
            || (type != elem.type)
            || (ghost_type != elem.ghost_type)
            || (kind != elem.kind));
  }

  bool operator<(const Element& rhs) const {
    bool res = (rhs == ElementNull) || ((this->kind < rhs.kind) ||
                                        ((this->kind == rhs.kind) &&
                                         ((this->ghost_type < rhs.ghost_type) ||
                                          ((this->ghost_type == rhs.ghost_type) &&
                                           ((this->type < rhs.type) ||
                                            ((this->type == rhs.type) &&
                                             (this->element < rhs.element)))))));
    return res;
  }

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
  
public:

  const ElementType & getType(){return type;}
  const UInt & getIndex(){return element;};
  const GhostType & getGhostType(){return ghost_type;}
  const ElementKind & getElementKind(){return kind;}

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

public:
  ElementType type;
  UInt element;
  GhostType ghost_type;
  ElementKind kind;
};

struct CompElementLess {
  bool operator() (const Element& lhs, const Element& rhs) const {
    return lhs < rhs;
  }
};

__END_AKANTU__

#endif /* __AKANTU_ELEMENT_HH__ */
