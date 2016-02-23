/*
    openDCM, dimensional constraint manager
    Copyright (C) 2015  Stefan Troeger <stefantroeger@gmx.net>

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

#include <boost/test/unit_test.hpp>
#include <boost/chrono.hpp>
#include <Eigen/Core>

#include "opendcm/core/kernel.hpp"
#include "opendcm/core/mapdatatype.hpp"

using namespace boost::chrono;

BOOST_AUTO_TEST_SUITE(Mapping_test_suit);

BOOST_AUTO_TEST_CASE(basic_math) {
    
       typedef dcm::Eigen3Kernel<double> Kernel;
       dcm::details::MapType<Kernel> t1, t2;
       double d1 = 2.;
       double d2 = 2.;
       
       //does assignement work?
       t1 = d1;
       t2 = d2;
       d1 = t2;
       
       //check if simple arithmetic operators work
       BOOST_CHECK_CLOSE(t1+t2, d1+d2, 1e-9);
       BOOST_CHECK_CLOSE(t1*t2, d1*d2, 1e-9);
       BOOST_CHECK_CLOSE(t1-t2, d1-d2, 1e-9);
       BOOST_CHECK_CLOSE(t1/t2, d1/d2, 1e-9);
       
       //are the io operators ok?
       //std::stringstream string;
       //string<<t1;
       //BOOST_CHECK_EQUAL(string.str(), std::string("2"));
       
       //check some simple arithmetics timings and see how they compare
#ifdef DCM_DEBUG
       int num = 1e6;
#else 
       int num = 1e8;
#endif
       volatile double t;
       system_clock::time_point time1 = system_clock::now();
       Eigen::VectorXd vec1(num), vec2(num);
       vec1.setRandom();
       vec2.setRandom();
       
       for(int i=0; i<num; ++i) {
           d1 = vec2(i);
           vec1(i) = d1*d2+d1+d2-(d1-d2)/(d1*d2);
       }
       
       system_clock::duration duration1 = system_clock::now() - time1;
      
       time1 = system_clock::now();
       for(int i=0; i<num; ++i) {
           t1 = vec2(i);
           vec1(i) = t1*t2+t1+t2-(t1-t2)/(t1*t2);
       }
       
       system_clock::duration duration2 = system_clock::now() - time1;       

       BOOST_CHECK_LT(duration_cast<microseconds>(duration2), 1.5*duration_cast<microseconds>(duration1));
      
       std::cout << "normal: " << duration_cast<milliseconds>(duration1)<<std::endl;
       std::cout << "custom: " << duration_cast<milliseconds>(duration2)<<std::endl;
}

BOOST_AUTO_TEST_CASE(eigen3_math) {
    
       typedef dcm::Eigen3Kernel<double> Kernel;
       typedef dcm::details::MapType<Kernel> MapType;
       
       Eigen::Matrix<MapType, 3, 1> t1(3., 2.,1.), t2(3.,2.,1.); //3 dimensional vectors
       Eigen::Matrix<double, 3, 1> d1(1.,2.,3.), d2(1.,2.,3.);
              
       //does assignement work?
       t1 = d1.cast<MapType>();
       t2 = d2.cast<MapType>();

      
       BOOST_CHECK_SMALL((t1.cast<double>() - d1).norm(), 1e-9);
       BOOST_CHECK_SMALL((t2.cast<double>() - d2).norm(), 1e-9);

                  
       //check if simple arithmetic operators work
       BOOST_CHECK_CLOSE(double((t1+t2).norm()), (d1+d2).norm(), 1e-9);
       BOOST_CHECK_CLOSE(double((t1*t2.transpose()).norm()), (d1*d2.transpose()).norm(), 1e-9);
       BOOST_CHECK_CLOSE(double((t1-t2).norm()), (d1-d2).norm(), 1e-9);
       BOOST_CHECK_CLOSE(double(t1.dot(t2)), d1.dot(d2), 1e-9);
       
       
       //check some simple arithmetics timings and see how they compare
#ifdef DCM_DEBUG
       int num = 1e5;
#else 
       int num = 1e8;
#endif
       volatile double t;
       Eigen::VectorXd vec1(num), vec2(num);
       vec1.setRandom();
       vec2.setRandom();
       
       system_clock::time_point time1 = system_clock::now();
       
       for(int i=0; i<num; ++i) {
           d1(1) = vec2(i);
           double tmp = (d1.dot(d2)*d1+d2-(d1-d2).norm()*(d1+3*d2)).norm();
           vec1(i) = tmp*(d1.cross(d2).norm() - (d1*d2.transpose()).norm());
       }
       
       system_clock::duration duration1 = system_clock::now() - time1;
      
       time1 = system_clock::now();
       for(int i=0; i<num; ++i) {
           t1(1) = vec2(i);
           double tmp = (t1.dot(t2)*t1+t2-(t1-t2).norm()*(t1+3*t2)).norm();
           vec1(i) = tmp*(t1.cross(t2).norm() - (t1*t2.transpose()).norm());
       }
       
       system_clock::duration duration2 = system_clock::now() - time1;       
   
       time1 = system_clock::now();
       for(int i=0; i<num; ++i) {
           t1(1) = vec2(i);
           double tmp = (t1.cast<double>().dot(d2)*t1.cast<double>()+d2-(t1.cast<double>()-d2).norm()*(t1.cast<double>()+3*d2)).norm();
           vec1(i) = tmp*(t1.cast<double>().cross(d2).norm() - (t1.cast<double>()*d2.transpose()).norm());
       }
       
       system_clock::duration duration3 = system_clock::now() - time1; 
       
       time1 = system_clock::now();
       for(int i=0; i<num; ++i) {
           t1(1) = vec2(i);
           double tmp = (t1.dot(d2.cast<MapType>())*t1+d2.cast<MapType>()-(t1-d2.cast<MapType>()).norm()*(t1+3*d2.cast<MapType>())).norm();
           vec1(i) = tmp*(t1.cross(d2.cast<MapType>()).norm() - (t1*d2.cast<MapType>().transpose()).norm());
       }
       
       system_clock::duration duration4 = system_clock::now() - time1;
       
       time1 = system_clock::now();
       Eigen::Vector3d v(1., 2., 3.);
       for(int i=0; i<(num-3); ++i) {
           d1(0) = vec2(i);
           d1(1) = vec2(i+1);
           d1(2) = vec2(i+2);
           vec1(i+1)= vec2(i);
           double tmp  = (d1.dot(d2)*d1+d2-(d1-d2).norm()*(d1+3*d2)).norm();
           vec1(i) = tmp*(d1.cross(d2).norm() - (d1*d2.transpose()).norm());
       }
       
       system_clock::duration duration5 = system_clock::now() - time1;
       
       time1 = system_clock::now();
       Eigen::Map<Eigen::Vector3d> m1(&d1(0)), m2(&d2(0));
       for(int i=0; i<num; ++i) {
           m1(1) = vec2(i);
           double tmp = (m1.dot(m2)*m1+m2-(m1-m2).norm()*(m1+3*m2)).norm();
           vec1(i) = tmp*(m1.cross(m2).norm() - (m1*m2.transpose()).norm());
       }
       
       system_clock::duration duration6 = system_clock::now() - time1;
       
       time1 = system_clock::now();
       for(int i=0; i<num; ++i) {
           m1(1) = vec2(i);
           double tmp = (m1.dot(d2)*m1+d2-(m1-d2).norm()*(m1+3*d2)).norm();
           vec1(i) = tmp*(m1.cross(d2).norm() - (m1*d2.transpose()).norm());
       }
       
       system_clock::duration duration7 = system_clock::now() - time1;
       
       std::cout << "normal:   " << duration_cast<milliseconds>(duration1)<<std::endl;
       std::cout << "custom:   " << duration_cast<milliseconds>(duration2)<<std::endl;
       std::cout << "mixed:    " << duration_cast<milliseconds>(duration3)<<std::endl;
       std::cout << "mixedI:   " << duration_cast<milliseconds>(duration4)<<std::endl;
       std::cout << "assigned: " << duration_cast<milliseconds>(duration5)<<std::endl;
       std::cout << "mapped:   " << duration_cast<milliseconds>(duration6)<<std::endl;
       std::cout << "mixedmap: " << duration_cast<milliseconds>(duration7)<<std::endl;
       
       BOOST_CHECK_LT(duration_cast<microseconds>(duration3), 1.5*duration_cast<microseconds>(duration1));
}

BOOST_AUTO_TEST_SUITE_END();
