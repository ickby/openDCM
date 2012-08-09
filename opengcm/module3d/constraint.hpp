/*
    openGCM, geometric constraint manager
    Copyright (C) 2012  Stefan Troeger <stefantroeger@gmx.net>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more detemplate tails.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef GCM_CONSTRAINT_H
#define GCM_CONSTRAINT_H

#include <boost/proto/core.hpp>
#include <boost/proto/domain.hpp>
#include <boost/proto/generate.hpp>

using namespace boost;

//allow all terminals and the addition of them
struct constraint_grammar
  : proto::or_<
        proto::plus< constraint_grammar, constraint_grammar >
      , proto::terminal< proto::_ >
    >
{};

//context to evaluate the terminals
struct constraint_context : proto::callable_context< constraint_context const >
{
    typedef void result_type;

    // Handle the placeholders:
    template<typename T>
    void operator()(proto::tag::terminal, single_constraint<T>) const
    {
        return this->args[I];
    }
};

// Forward-declare an expression wrapper
template<typename Expr>
struct constraint;

// Define a constraint domain. Expression within
// the constraint domain will be wrapped in the
// constraint<> expression wrapper.
struct constraint_domain
  : proto::domain< proto::generator< constraint<> > >
{};

// Define a constraint expression wrapper. It behaves just like
// the expression it wraps, but with an extra operator() member
// function that evaluates the expression.    
template<typename Expr>
struct constraint
  : proto::extends<Expr, constraint<Expr>, constraint_domain>
{
    typedef
        proto::extends<Expr, constraint<Expr>, constraint_domain>
    base_type;

    constraint(Expr const &expr = Expr())
      : base_type(expr)
    {}

    typedef double result_type;

    // Overload operator() to invoke proto::eval() with
    // our constraint_context.
    template<typename Vector>
    double operator()(Vector v) const
    {
        constraint_context ctx;
        ctx.args.push_back(a1);
        ctx.args.push_back(a2);

        return proto::eval(*this, ctx);
    }
};

#endif