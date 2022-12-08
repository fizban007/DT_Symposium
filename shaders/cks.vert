/*
 * Copyright (C) 2016 Chi-kwan Chan
 * Copyright (C) 2016 Steward Observatory
 *
 * This file is part of GRay2.
 *
 * GRay2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GRay2 is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GRay2.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \file
 ** Cartesian Kerr-Schild coordinate specific schemes
 **
 ** Implement the coordinate specific functions getuu(), icond(), and
 ** rhs() in the Cartesian form of the Kerr-Schild coordiantes.  Let
 ** \f$t\f$, \f$x\f$, \f$y\f$, \f$z\f$ be the coordinates, the
 ** Cartesian Kerr-Schild metric is given by
 ** \f[
 **   g_{\mu\nu} = \gamma_{\mu\nu} + f l_\mu l_\nu
 ** \f]
 ** where \f$\gamma_{\mu\nu}\f$ is the Minkowski metric, \f$f\f$ and
 ** \f$l_\mu\f$ are defined by
 ** \f[
 **   f = \frac{2r^3}{r^4 + a^2 z^2} \mbox{ and }
 **   l_\mu = \left(1, \frac{rx + ay}{r^2 + a^2},
 **                    \frac{ry - ax}{r^2 + a^2},
 **                    \frac{z}{r}\right),
 ** \f]
 ** respectively, and \f$r\f$ is defined implicitly by\f$ x^2 + y^2 +
 ** z^2 = r^2 + a^2 (1 - z^2 / r^2)\f$.
 **/

/**
 ** The sqaure of vector u at the spacetime event q in Kerr-Schild coordiantes
 **
 ** Compute \f$u\cdot u \equiv g_{\alpha\beta} u^\alpha u^\beta\f$,
 ** where \f$g_{\alpha\beta}\f$ is the Cartesian form of the
 ** Kerr-Schild metric.
 **
 ** \return The square of u at q
 **/
static double
getuu(double4 q, /**< Spacetime event "location" */
      double4 u) /**< The vector being squared   */
{
  double  aa = a_spin * a_spin;
  double  zz = q.s3 * q.s3;
  double  kk = 0.5 * (q.s1 * q.s1 + q.s2 * q.s2 + zz - aa);
  double  dd = sqrt(kk * kk + aa * zz);
  double  rr = dd + kk;
  double  r  = sqrt(rr);

  double  f  = 2.0 * rr * r / (rr * rr + aa * zz);
  double  lx = (r * q.s1 + a_spin * q.s2) / (rr + aa);
  double  ly = (r * q.s2 - a_spin * q.s1) / (rr + aa);
  double  lz = q.s3 / r;

  double4 gt = {-1. + f*1.*1.,      f*1.*lx,      f*1.*ly,      f*1.*lz};
  double4 gx = {      f*lx*1., 1. + f*lx*lx,      f*lx*ly,      f*lx*lz};
  double4 gy = {      f*ly*1.,      f*ly*lx, 1. + f*ly*ly,      f*ly*lz};
  double4 gz = {      f*lz*1.,      f*lz*lx,      f*lz*ly, 1. + f*lz*lz};

  return (dot(gt, u) * u.s0 +
          dot(gx, u) * u.s1 +
          dot(gy, u) * u.s2 +
          dot(gz, u) * u.s3);
}

/**
 ** The initial conditions of a ray in an image plane
 **
 ** To perform ray tracing calculations of an image in Kerr spacetime,
 ** we follow Johannsen & Psaltis (2010) and consider an observer
 ** viewing the central black hole from a large distance \p r_obs and
 ** at an inclination angle \p i_obs from its rotation axis (see
 ** Figure 1 of Psaltis & Johannsen 2012).  We set up a virtual image
 ** plane that is perpendicular to the line of sight and centered at
 ** \f$\phi\f$ = \p j_obs of the spacetime.  We define the set of
 ** local Cartesian coordinates (\p alpha, \p beta) on the image plane
 ** such that the \p beta axis is along the same fiducial plane and
 ** the \p alpha axis is perpendicular to it.  These input parameters
 ** define a unique ray, whose initial spacetime position and
 ** wavevector are then computed by icond().
 **
 ** \todo Consider the fact that Kerr spacetime is not flat
 **
 ** \return The initial conditions of a ray
 **/
static double8
icond(double r_obs, /**< Distance of the observer from the black hole */
      double i_obs, /**< Inclination angle of the observer in degrees */
      double j_obs, /**< Azimuthal   angle of the observer in degrees */
      double alpha, /**< One of the local Cartesian coordinates       */
      double beta)  /**< The other  local Cartesian coordinate        */
{
  double  ci, si = sincos(M_PI * i_obs / 180.0, &ci);
  double  cj, sj = sincos(M_PI * j_obs / 180.0, &cj);

  double  R = r_obs * si - beta  * ci; /* cylindrical radius */
  double  z = r_obs * ci + beta  * si;
  double  y = R     * sj - alpha * ci;
  double  x = R     * cj + alpha * sj;

  double4 q = (double4){0.0, x, y, z};
  double4 u = (double4){1.0, si * cj, si * sj, ci};

  return (double8){q, u};
}

/**
 ** The right hand sides of the geodesic equations in Kerr-Schild coordiantes
 **
 ** One of the breakthroughs we achieve in GRay2 is that, by a series
 ** of mathematical manipulations and regrouping, we significantly
 ** reduce the operation count of the geodesic equations in the
 ** Cartesian Kerr-Schild coordinates.  Let \f$\lambda\f$ be the
 ** affine parameter and \f$\dot{x}^\mu \equiv dx^\mu/d\lambda\f$.  We
 ** show in Chan et al. (2017) that the geodesic equations in the
 ** Cartesian KS coordinates can be optimized to the following form:
 ** \f[
 **  \ddot{x}^\mu = - \left(\eta^{\mu\beta} \dot{x}^\alpha -
 **                         \frac{1}{2}\eta^{\mu\alpha} \dot{x}^\beta\right)
 **                 \dot{x}_{\beta,\alpha} + F l^\mu
 ** \f]
 ** where
 ** \f[
 **   F = f \left(l^\beta \dot{x}^\alpha -
 **               \frac{1}{2}l^\alpha \dot{x}^\beta\right)
 **       \dot{x}_{\beta,\alpha}.
 ** \f]
 ** In this new form, the right hand sides (RHS) of the geodesic
 ** equations have only ~65% more floating-point operations than in
 ** the Boyer-Lindquist coordinates.  Furthermore, the evaluation of
 ** the RHS uses many matrix-vector products, which are optimized in
 ** modern hardwares.
 **
 ** \return The right hand sides of the geodesic equations
 **/
static double8
rhs(double8 s) /**< State of the ray */
{
  double4 q = s.s0123;
  double4 u = s.s4567;

  double f,  dx_f,  dy_f,  dz_f;
  double lx, dx_lx, dy_lx, dz_lx;
  double ly, dx_ly, dy_ly, dz_ly;
  double lz, dx_lz, dy_lz, dz_lz;

  double  hDxu, hDyu, hDzu;
  double4 uD;
  double  tmp;

  {
    double dx_r, dy_r, dz_r;
    double r, ir, iss;
    {
      double aa = a_spin * a_spin;
      double rr;
      {
        double zz = q.s3 * q.s3;
        double dd;
        {
          double kk = 0.5 * (q.s1 * q.s1 + q.s2 * q.s2 + zz - aa);
          dd = sqrt(kk * kk + aa * zz);
          rr = dd + kk;
        }
        r  = sqrt(rr);
        ir = 1.0 / r;
        {
          double ss = rr + aa;
          iss  = 1.0 / ss;
          tmp  = 0.5 / (r * dd);
          dz_r = tmp * ss * q.s3;
          tmp *= rr;
        }
        dy_r = tmp * q.s2;
        dx_r = tmp * q.s1;
        tmp  = 2.0 / (rr + aa * zz / rr);
      }
      f    = tmp *  r;
      dx_f = tmp *  dx_r * (3.0 - 2.0 * rr * tmp);
      dy_f = tmp *  dy_r * (3.0 - 2.0 * rr * tmp);
      dz_f = tmp * (dz_r * (3.0 - 2.0 * rr * tmp) - tmp * aa * q.s3 * ir);
    } /* 48 (-8) FLOPs; estimated FLoating-point OPerations, the number
         in the parentheses is (the negative of) the number of FMA */
    {
      double m2r  = -2.0 * r;
      double issr =  iss * r;
      double issa =  iss * a_spin;

      lx    = iss * (q.s1 * r + q.s2 * a_spin);
      tmp   = iss * (q.s1 + m2r * lx);
      dx_lx = tmp * dx_r + issr;
      dy_lx = tmp * dy_r + issa;
      dz_lx = tmp * dz_r;

      ly    = iss * (q.s2 * r - q.s1 * a_spin);
      tmp   = iss * (q.s2 + m2r * ly);
      dx_ly = tmp * dx_r - issa;
      dy_ly = tmp * dy_r + issr;
      dz_ly = tmp * dz_r;

      lz    = q.s3 * ir;
      tmp   = -lz * ir;
      dx_lz = tmp * dx_r;
      dy_lz = tmp * dy_r;
      dz_lz = tmp * dz_r + ir;
    } /* 35 (-9) FLOPs */
  }

  {
    double  flu;
    double4 Dx, Dy, Dz;
    {
      double lu = u.s0 + lx * u.s1 + ly * u.s2 + lz * u.s3;
      flu   = f * lu;
      Dx.s0 = dx_f * lu + f * (dx_lx * u.s1 + dx_ly * u.s2 + dx_lz * u.s3);
      Dy.s0 = dy_f * lu + f * (dy_lx * u.s1 + dy_ly * u.s2 + dy_lz * u.s3);
      Dz.s0 = dz_f * lu + f * (dz_lx * u.s1 + dz_ly * u.s2 + dz_lz * u.s3); /* 31 (-12) FLOPs */
    }
    Dx.s1 = Dx.s0 * lx + flu * dx_lx;
    Dx.s2 = Dx.s0 * ly + flu * dx_ly;
    Dx.s3 = Dx.s0 * lz + flu * dx_lz; /* 9 (-3) FLOPs */

    Dy.s1 = Dy.s0 * lx + flu * dy_lx;
    Dy.s2 = Dy.s0 * ly + flu * dy_ly;
    Dy.s3 = Dy.s0 * lz + flu * dy_lz; /* 9 (-3) FLOPs */

    Dz.s1 = Dz.s0 * lx + flu * dz_lx;
    Dz.s2 = Dz.s0 * ly + flu * dz_ly;
    Dz.s3 = Dz.s0 * lz + flu * dz_lz; /* 9 (-3) FLOPs */

    hDxu = 0.5 * dot(Dx, u);
    hDyu = 0.5 * dot(Dy, u);
    hDzu = 0.5 * dot(Dz, u); /* 24 (-9) FLOPs */

    uD  = u.s1 * Dx + u.s2 * Dy + u.s3 * Dz; /* 20 (-8) FLOPs */

    tmp = f * (-uD.s0 + lx * (uD.s1 - hDxu) + ly * (uD.s2 - hDyu) + lz * (uD.s3 - hDzu)); /* 10 (-3) FLOPs */
  }

  return (double8){u,
        uD.s0 -      tmp,
        hDxu - uD.s1 + lx * tmp,
        hDyu - uD.s2 + ly * tmp,
        hDzu - uD.s3 + lz * tmp}; /* 10 (-3) FLOPs */
}
