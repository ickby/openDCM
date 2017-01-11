/*
    openDCM, dimensional constraint manager
    Copyright (C) 2017  Stefan Troeger <stefantroeger@gmx.net>

    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License along
    with this library; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

/** @addtogroup Core
 * @{*/

/** @defgroup Reduction Constraint system reduction
 *
 * @brief Provides ways to recombine symbolic geometries and constraints
 * 
 * \section reductiongraph Purpose and scope
 * 
 * A geometric constraint system as setup by the user is an intutively created network of geometries 
 * and constraint, but it is not computational optimal. Many of the individual created combinations 
 * can be combined into annother setup with less equations and less free parameters. The redction 
 * infrastructure provided by dcm is the way to achieve this recombination and to create a computational
 * optimal representation of the system. 
 *  
 * \subsection structure The reduction graph structure
 * If a geometry is connected with constraints to annother one it needs to be decided if it can be 
 * represented as a function of the base geometry. 
 * \image html reductiongraph.svg A node / edge reduction structure
 * 
 * 
 *@}
 */

