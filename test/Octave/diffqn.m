
function [m_rot, m_diffrot] = diffqn(m_normQ)

	nor = norm(m_normQ);
        Q = [m_normQ(1)/nor, m_normQ(2)/nor, m_normQ(3)/nor, nor];
	Q = Q./norm(Q);

 %       /* now calculate the gradient quaternions and calculate the diff rotation matrices
 %        * normQ = (a,b,c)
 %        * sn0 = ||normQ||^2, n0 = ||normQ||
 %        *
 %        * Q = (a/n0, b/n0, c/n0, n0)
 %        * ||Q|| = n = 1/n0 * sqrt( a^2+b^2+c^2+sn0^2 )
 %        * n2 = sqrt( a^2+b^2+c^2+sn0^2 )
 %        *
 %        * unit Quaternion uQ = (x y z w) = (a/(n0*n), b/(n0*n), c/(n0*n), n0/n)
 %        * uQ = (a/n2, b/n2, c/n2, sn0/n2)
 %        */

        sn0 = norm(m_normQ)^2;
        n0 = norm(m_normQ);
        n2 = sqrt(sn0 + sn0^2);

        %// d(1/n2)/dx and dy and dz
        ddn2 = -(1+2*sn0)/(n2^3);
        n2_dda = m_normQ(1)*ddn2;
        n2_ddb = m_normQ(2)*ddn2;
        n2_ddc = m_normQ(3)*ddn2;

        %//dxa = da/dx
        dxa = 1/n2 + m_normQ(1)*n2_dda;
        dxb = m_normQ(1)*n2_ddb;
        dxc = m_normQ(1)*n2_ddc;

        dya = m_normQ(2)*n2_dda;
        dyb = 1/n2 + m_normQ(2)*n2_ddb;
        dyc = m_normQ(2)*n2_ddc;

        dza = m_normQ(3)*n2_dda;
        dzb = m_normQ(3)*n2_ddb;
        dzc = 1/n2 + m_normQ(3)*n2_ddc;

        dwa = 2*m_normQ(1)/n2 + sn0*n2_dda;
        dwb = 2*m_normQ(2)/n2 + sn0*n2_ddb;
        dwc = 2*m_normQ(3)/n2 + sn0*n2_ddc;
	
	%//write in the diffrot matrix, starting with duQ/dx
        m_diffrot(1,1) = -4.0*(Q(2)*dya+Q(3)*dza);
        m_diffrot(1,2) = -2.0*(Q(4)*dza+dwa*Q(3))+2.0*(Q(1)*dya+dxa*Q(2));
        m_diffrot(1,3) = 2.0*(dwa*Q(2)+Q(4)*dya)+2.0*(dxa*Q(3)+Q(1)*dza);
        m_diffrot(2,1) = 2.0*(Q(4)*dza+dwa*Q(3))+2.0*(Q(1)*dya+dxa*Q(2));
        m_diffrot(2,2) = -4.0*(Q(1)*dxa+Q(3)*dza);
        m_diffrot(2,3) = -2.0*(dwa*Q(1)+Q(4)*dxa)+2.0*(dya*Q(3)+Q(2)*dza);
        m_diffrot(3,1) = -2.0*(dwa*Q(2)+Q(4)*dya)+2.0*(dxa*Q(3)+Q(1)*dza);
        m_diffrot(3,2) = 2.0*(dwa*Q(1)+Q(4)*dxa)+2.0*(dya*Q(3)+Q(2)*dza);
        m_diffrot(3,3) = -4.0*(Q(1)*dxa+Q(2)*dya);

        m_diffrot(1,4) = -4.0*(Q(2)*dyb+Q(3)*dzb);
        m_diffrot(1,5) = -2.0*(Q(4)*dzb+dwb*Q(3))+2.0*(Q(1)*dyb+dxb*Q(2));
        m_diffrot(1,6) = 2.0*(dwb*Q(2)+Q(4)*dyb)+2.0*(dxb*Q(3)+Q(1)*dzb);
        m_diffrot(2,4) = 2.0*(Q(4)*dzb+dwb*Q(3))+2.0*(Q(1)*dyb+dxb*Q(2));
        m_diffrot(2,5) = -4.0*(Q(1)*dxb+Q(3)*dzb);
        m_diffrot(2,6) = -2.0*(dwb*Q(1)+Q(4)*dxb)+2.0*(dyb*Q(3)+Q(2)*dzb);
        m_diffrot(3,4) = -2.0*(dwb*Q(2)+Q(4)*dyb)+2.0*(dxb*Q(3)+Q(1)*dzb);
        m_diffrot(3,5) = 2.0*(dwb*Q(1)+Q(4)*dxb)+2.0*(dyb*Q(3)+Q(2)*dzb);
        m_diffrot(3,6) = -4.0*(Q(1)*dxb+Q(2)*dyb);

        m_diffrot(1,7) = -4.0*(Q(2)*dyc+Q(3)*dzc);
        m_diffrot(1,8) = -2.0*(Q(4)*dzc+dwc*Q(3))+2.0*(Q(1)*dyc+dxc*Q(2));
        m_diffrot(1,9) = 2.0*(dwc*Q(2)+Q(4)*dyc)+2.0*(dxc*Q(3)+Q(1)*dzc);
        m_diffrot(2,7) = 2.0*(Q(4)*dzc+dwc*Q(3))+2.0*(Q(1)*dyc+dxc*Q(2));
        m_diffrot(2,8) = -4.0*(Q(1)*dxc+Q(3)*dzc);
        m_diffrot(2,9) = -2.0*(dwc*Q(1)+Q(4)*dxc)+2.0*(dyc*Q(3)+Q(2)*dzc);
        m_diffrot(3,7) = -2.0*(dwc*Q(2)+Q(4)*dyc)+2.0*(dxc*Q(3)+Q(1)*dzc);
        m_diffrot(3,8) = 2.0*(dwc*Q(1)+Q(4)*dxc)+2.0*(dyc*Q(3)+Q(2)*dzc);
        m_diffrot(3,9) = -4.0*(Q(1)*dxc+Q(2)*dyc);


  tx  = 2*Q(1);
  ty  = 2*Q(2);
  tz  = 2*Q(3);
  twx = tx*Q(4);
  twy = ty*Q(4);
  twz = tz*Q(4);
  txx = tx*Q(1);
  txy = ty*Q(1);
  txz = tz*Q(1);
  tyy = ty*Q(2);
  tyz = tz*Q(2);
  tzz = tz*Q(3);

  m_rot(1,1) = 1-(tyy+tzz);
  m_rot(1,2) = txy-twz;
  m_rot(1,3) = txz+twy;
  m_rot(2,1) = txy+twz;
  m_rot(2,2) = 1-(txx+tzz);
  m_rot(2,3) = tyz-twx;
  m_rot(3,1) = txz-twy;
  m_rot(3,2) = tyz+twx;
  m_rot(3,3) = 1-(txx+tyy);
 
endfunction