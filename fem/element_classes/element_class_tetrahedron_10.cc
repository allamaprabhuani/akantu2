/**
 * @file   element_class_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 10:28:28 2010
 *
 * @brief  Specialization of the element_class class for the type _tetrahedron_10
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 * @section DESCRIPTION
 *
 * @verbatim
	\zeta
	  ^
	  |
       (0,0,1)
	  x
	  |` .
	  |  `  .
	  |    `   .
	  |      `    .  (0,0.5,0.5)
	  |        `    x.
	  |     q4 o `      .                   \eta
	  |            `       .             -,
(0,0,0.5) x             ` x (0.5,0,0.5)  -
	  |                `        x-(0,1,0)
	  |              q3 o`   -   '
	  |        0,0.5,0)  - `      '
	  |             x-       `     x (0.5,0.5,0)
	  |     q1 o -         o q2`    '
	  |      -                   `   '
	  |  -                         `  '
	  x---------------x--------------` x-----> \xi
       (0,0,0)        (0.5,0,0)        (1,0,0)
 @endverbatim
 *
 * @subsection coords Nodes coordinates
 *
 * @f[
 * \begin{array}{lll}
 *   \xi_{0}  = 0   &  \eta_{0}  = 0   &  \zeta_{0}  = 0   \\
 *   \xi_{1}  = 1   &  \eta_{1}  = 0   &  \zeta_{1}  = 0   \\
 *   \xi_{2}  = 0   &  \eta_{2}  = 1   &  \zeta_{2}  = 0   \\
 *   \xi_{3}  = 0   &  \eta_{3}  = 0   &  \zeta_{3}  = 1   \\
 *   \xi_{4}  = 1/2 &  \eta_{4}  = 0   &  \zeta_{4}  = 0   \\
 *   \xi_{5}  = 1/2 &  \eta_{5}  = 1/2 &  \zeta_{5}  = 0   \\
 *   \xi_{6}  = 0   &  \eta_{6}  = 1/2 &  \zeta_{6}  = 0   \\
 *   \xi_{7}  = 0   &  \eta_{7}  = 0   &  \zeta_{7}  = 1/2 \\
 *   \xi_{8}  = 1/2 &  \eta_{8}  = 0   &  \zeta_{8}  = 1/2 \\
 *   \xi_{9}  = 0   &  \eta_{9}  = 1/2 &  \zeta_{9}  = 1/2
 * \end{array}
 * @f]
 *
 * @subsection shapes Shape functions
 * @f[
 * \begin{array}{llll}
 *     N1  = (1 - \xi - \eta - \zeta) (1 - 2 \xi - 2 \eta - 2 \zeta)
 *           & \frac{\partial N1}{\partial \xi}    = 
 *           & \frac{\partial N1}{\partial \eta}   =
 *           & \frac{\partial N1}{\partial \zeta}  = \\
 *     N2  = \xi (2 \xi - 1)
 *           & \frac{\partial N2}{\partial \xi}    =
 *           & \frac{\partial N2}{\partial \eta}   = 0
 *           & \frac{\partial N2}{\partial \zeta}  = 0 \\
 *     N3  = \eta (2 \eta - 1)
 *           & \frac{\partial N3}{\partial \xi}    = 0
 *           & \frac{\partial N3}{\partial \eta}   =
 *           & \frac{\partial N3}{\partial \zeta}  = 0 \\
 *     N4  = \zeta (2 \zeta - 1)
 *           & \frac{\partial N4}{\partial \xi}    = 0
 *           & \frac{\partial N4}{\partial \eta}   = 0
 *           & \frac{\partial N4}{\partial \zeta}  =  \\
 *     N5  = 4 \xi (1 - \xi - \eta - \zeta)
 *           & \frac{\partial N5}{\partial \xi}    = 
 *           & \frac{\partial N5}{\partial \eta}   = -4 \xi
 *           & \frac{\partial N5}{\partial \zeta}  = -4 \xi \\
 *     N6  = 4 \xi \eta
 *           & \frac{\partial N6}{\partial \xi}    = 4 \eta
 *           & \frac{\partial N6}{\partial \eta}   = 4 \xi
 *           & \frac{\partial N6}{\partial \zeta}  = 0 \\
 *     N7  = 4 \eta (1 - \xi - \eta - \zeta)
 *           & \frac{\partial N7}{\partial \xi}    = -4 \xi
 *           & \frac{\partial N7}{\partial \eta}   =
 *           & \frac{\partial N7}{\partial \zeta}  = -4 \xi \\
 *     N8  = 4 \zeta (1 - \xi - \eta - \zeta)
 *           & \frac{\partial N8}{\partial \xi}    = -4 \zeta
 *           & \frac{\partial N8}{\partial \eta}   = -4 \zeta
 *           & \frac{\partial N8}{\partial \zeta}  = \\
 *     N9  = 4 \zeta \xi
 *           & \frac{\partial N9}{\partial \xi}    = 4 \zeta
 *           & \frac{\partial N9}{\partial \eta}   = 0
 *           & \frac{\partial N9}{\partial \zeta}  = 4 \xi \\
 *     N10 = 4 \eta \zeta
 *           & \frac{\partial N10}{\partial \xi}   = 0
 *           & \frac{\partial N10}{\partial \eta}  = 4 \zeta
 *           & \frac{\partial N10}{\partial \zeta} = 4 \eta \\
 * \end{array}
 * @f]
 *
 * @subsection quad_points Position of quadrature points
 * @f[
 * a = \frac{5 - \sqrt{5}}{20}\\
 * b = \frac{5 + 3 \sqrt{5}}{20}
 * \begin{array}{lll}
 *   \xi_{q_0}  = a   &  \eta_{q_0}  = a   &  \zeta_{q_0}  = a \\
 *   \xi_{q_1}  = b   &  \eta_{q_1}  = a   &  \zeta_{q_1}  = a \\
 *   \xi_{q_2}  = a   &  \eta_{q_2}  = b   &  \zeta_{q_2}  = a \\
 *   \xi_{q_3}  = a   &  \eta_{q_3}  = a   &  \zeta_{q_3}  = b
 * \end{array}
 * @f]
 */

/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
