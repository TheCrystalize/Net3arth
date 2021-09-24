/* helper functions */

function C(real, imaginary) {
  return {
    re: real,
    im: imaginary
  };
}

function conj(z) {
  return {
    re: z.re,
    im: -z.im
  }
}

function div(z, c) {
  const s = 1 / (c.re * c.re + c.im * c.im);
  return {
    re: (z.re * c.re + z.im * c.im) * s,
    im: (z.im * c.re - z.re * c.im) * s
  }
}

function divScalar(z, n) {
  return {
    re: z.re / n,
    im: z.im / n
  }
}

function add(z, c) {
  return {
    re: z.re + c.re,
    im: z.im + c.im
  }
}

function neg(z) {
  return {
    re: -z.re,
    im: -z.im
  }
}

function sub(z, c) {
  return {
    re: z.re - c.re,
    im: z.im - c.im
  }
}

function addScalar(z, s) {
  return {
    re: z.re + s,
    im: z.im
  }
}

function mult(z, c) {
  return {
    re: z.re * c.re - z.im * c.im,
    im: z.re * c.im + z.im * c.re
  }
}

function multScalar(z, s) {
  return {
    re: z.re * s,
    im: z.im * s
  }
}

function sqrt(z) {
  const s = Math.hypot(z.re, z.im),
    sgn = z.im < 0 ? -1 : 1;
  return multScalar({
    re: Math.sqrt(s + z.re),
    im: sgn * Math.sqrt(s - z.re)
  }, 0.5 * Math.SQRT2);
}

function log(z) {
  return {
    re: 0.5 * Math.log(z.re * z.re + z.im * z.im),
    im: Math.atan2(z.im, z.re)
  }
}

function pow(z, p) {
  const n = p * Math.atan2(z.im, z.re);
  return {
    re: Math.cos(n),
    im: Math.sin(n) * Math.exp(p * Math.log(z.re * z.re + z.im * z.im))
  }
}

function exp(z) {
  const e = Math.exp(z.re);
  return {
    re: Math.cos(z.im) * e,
    im: Math.sin(z.im) * e
  }
}

function dot(z, c) {
  return (z.re * c.re + z.im * c.im)
}

function sinh(z) {
  return add(divScalar(exp(z), 2), neg(divScalar(exp(neg(z)), 2)))
}

function cosh(z) {
  return add(divScalar(exp(neg(z)), 2), divScalar(exp(z), 2))
}

function tanh(z) {
  return div(sinh(z), cosh(z))
}

function gaussRnd() {
  return (Math.random() + Math.random()) * 2.0 - 2.0;
}

/* transforms */

function reset() {
  return z=>{
    return {
      re: 0.000001,
      im: 0.000002,
      red: 1,
      green: 1,
      blue: 1,
      alpha: 1,
    }
  }
}

function arcsinh() {
  return z => {
    return {
      ...z,
      ...multScalar(
        log(
          add(
            z,
            sqrt(addScalar(mult(z, z), 1))
          )
        ),
        2 / Math.PI)
    }
  }
}

function arctanh() {
  return z => {
    return {
      ...z,
      ...multScalar(
        log(
          div(
            addScalar(z, 1),
            addScalar(neg(z), 1)
          )
        ),
        1 / Math.PI)
    }
  }
}

function blurCircle() {
  return z => {
    let a = Math.random() * Math.PI * 2,
      r = Math.sqrt(Math.random());
    return {
      ...z,
      re: Math.cos(a) * r,
      im: Math.sin(a) * r
    }
  }
}

function blurGasket() {
  return z => {
    let a = Math.random() * Math.PI * 2,
      r = 1 / Math.sqrt(Math.abs(Math.random() - Math.random()));
    return {
      ...z,
      re: Math.random() - 0.5,
      im: Math.sin(a) * r
    }
  }
}

function blurGaussian(pow) {
  return z => {
    let a = gaussRnd() * Math.PI * 2,
      s = Math.sqrt(Math.abs(gaussRnd())) * pow;
    return {
      ...z,
      re: Math.cos(a) * s + z.re,
      im: Math.sin(a) * s + z.im
    }
  }
}

function blurSine(pow) {
  return z => {
    let a = Math.random() * 2 * Math.PI,
      u = Math.random();
    let r = (pow == 1 ? Math.acos(u * 2 - 1) : Math.acos(Math.exp(Math.log(1 - u) * power) * 2 - 1)) / Math.PI;
    return {
      ...z,
      re: Math.cos(a) * r,
      im: Math.sin(a) * r
    }
  }
}

function blurSquare() {
  return z => {
    return {
      ...z,
      re: Math.random() - 0.5,
      im: Math.random() - 0.5
    }
  }
}

function bubble() {
  return z => {
    let r = (dot(z, z) + 4)
    return {
      ...z,
      re: z.re / r,
      im: -z.im / r
    }
  }
}

function circleInv() {
  return z => {
    let r = 1 / dot(z, z)
    return {
      ...z,
      re: z.re * r,
      im: -z.im * r
    }
  }
}

function cpow(rotation, depth, divisor, spread) {
  let ang = 2 * Math.PI / divisor,
    a1 = depth < 0 ? -Math.log(-depth) * rotation : Math.log(depth) * rotation;
  let a = Math.atan2(a1, 2 * Math.PI);
  let c = Math.cos(a) * rotation * Math.cos(a) / divisor,
    d = Math.cos(a) * rotation * Math.sin(a) / divisor;
  let coeff = d == 0 ? 0 : -0.095 * spread / d;
  return z => {
    let a0 = Math.atan2(z.im, z.re);
    let a1 = a0 < 0 ? a0 + 2 * Math.PI : a0;
    let a2 = Math.cos(a1 * 0.5) < Math.random() * 2 - 1 ? a1 - 2 * Math.PI : a1;
    let a = a2 + ((Math.random() > 0.5 ? 2 : -2) * Math.PI) * Math.floor(Math.log(Math.random()) * coeff + 0.5),
      lnr2 = Math.log(dot(z, z));
    let r = Math.exp(c * 0.5 * lnr2 - d * a),
      angle = c * a + d * 0.5 * lnr2 + ang * Math.floor(Math.random() * 32767);
    return {
      ...z,
      re: Math.sin(angle) * r,
      im: Math.cos(angle) * r
    }
  }
}

function cylinder() {
  return z => {
    return {
      ...z,
      re: z.re / Math.sqrt(z.re * z.re + 1),
      im: z.im
    }
  }
}

function disc() {
  return z => {
    let r = Math.sqrt(dot(z, z)) * Math.PI,
      a = Math.atan2(z.im, z.re) / Math.PI;
    return {
      ...z,
      re: Math.sin(r) * a,
      im: Math.cos(r) * a
    }
  }
}

function hypershape(n) {
  let alpha = 2 * Math.PI / n;
  let alphacoeff = Math.tan(alpha * 0.5) * 2,
    beta = Math.SQRT2 * Math.cos(alpha * 0.5);
  return z => {
    let da = (Math.atan2(z.im, z.re) + Math.PI) / alpha,
      rad = Math.sqrt(dot(z, z));
    let za0 = Math.floor(da);
    let xa0 = da - za0,
      xa1, za, si;
    if (xa0 > 0.5) {
      xa1 = 1 - xa0;
      za = za0 + 1;
      si = -1;
    } else {
      xa1 = xa0;
      za = za0;
      si = 1;
    }
    let xa = Math.atan(xa1 * alphacoeff) / alpha;
    let co = beta / Math.cos(xa * alpha),
      ang = (za + si * xa) * alpha - Math.PI;

    let p = multScalar({
      re: Math.cos(ang),
      im: Math.sin(ang)
    }, co * rad);
    let r = dot(p, p);

    return {
      ...z,
      ...(divScalar(multScalar(p, Math.SQRT2 - Math.sqrt(2 - r)), r))
    }
  }
}

function hypershift(p) {
  return z => {
    return {
      ...z,
      ...div(add(z, p), addScalar(mult(conj(p), z), 1))
    }
  }
}

function hypertile3(p, q, r, shift) {
  let rad = Math.PI / 180;

  let o = Math.acosh((Math.cos(Math.PI / p) + Math.cos(Math.PI / q) * Math.cos(Math.PI / r)) / (Math.sin(Math.PI / q) * Math.sin(Math.PI / r)));
  let a = Math.asinh(Math.sin(Math.PI / q) / Math.sin(Math.PI / p) * Math.sinh(o));
  let b = Math.asinh(Math.sin(Math.PI / r) * Math.sinh(a)),
    c = Math.asinh(Math.sin(Math.PI / r) * Math.sinh(a) / Math.sin(Math.PI / q));
  let h = Math.tanh(b / 2),
    b1 = Math.tanh(Math.acosh(Math.cosh(c) / Math.cosh(b)) / 2),
    b2 = Math.tanh(Math.acosh(Math.cosh(a) / Math.cosh(b)) / 2),
    rot1 = 360 / p * rad,
    rot2 = 360 / q * rad,
    rot3 = 360 / r * rad;
  let r01 = C(Math.cos(rot1), Math.sin(rot1)),
    r02 = C(Math.cos(rot2), Math.sin(rot2)),
    r03 = C(Math.cos(rot3), Math.sin(rot3));

  const c1 = C(h, 0),
    c2 = C(0, b1),
    c3 = C(0, b2);

  return z => {
    let n, pfr, z0;
    if (shift < 0.25) {
      n = 0;
      pfr = 0;
      z0 = z;
    } else if (shift < 0.5) {
      n = p;
      pfr = rot1;
      z0 = div(add(z, c1), addScalar(mult(c1, z), 1));
    } else if (shift < 0.75) {
      n = q;
      pfr = rot2;
      z0 = div(add(z, neg(c2)), addScalar(mult(c2, z), 1));
    } else {
      n = r;
      pfr = rot3;
      z0 = div(add(z, c3), addScalar(mult(neg(c3), z), 1));
    }

    let m0 = div(add(z0, neg(c1)), addScalar(mult(neg(c1), z0), 1));
    let r1 = mult(r01, m0);
    let m0f = div(add(r1, c1), addScalar(mult(c1, r1), 1)),
      m1 = div(add(z0, c2), addScalar(mult(neg(c2), z0), 1));
    let r2 = mult(r02, m1);
    let m1f = div(add(r2, neg(c2)), addScalar(mult(c2, r2), 1)),

      m2 = div(add(z0, neg(c3)), addScalar(mult(c3, z0), 1));
    let r3 = mult(r03, m2);
    let m2f = div(add(r3, c3), addScalar(mult(neg(c3), r3), 1));

    let fr = Math.floor(Math.random() * n) * pfr,
      rnd = Math.random(),
      f3, f0, f;

    if (rnd < 1 / 3) {
      f3 = m0f
    } else if (rnd < 2 / 3) {
      f3 = m1f
    } else {
      f3 = m2f
    }

    if (shift < 0.25) {
      f0 = f3
    } else if (shift < 0.5) {
      f0 = div(add(f3, neg(c1)), addScalar(mult(neg(c1), f3), 1))
    } else if (shift < 0.75) {
      f0 = div(add(f3, c2), addScalar(mult(neg(c2), f3), 1))
    } else {
      f0 = div(add(f3, neg(c3)), addScalar(mult(c3, f3), 1))
    }

    if (shift < 0.25) {
      f = f0
    } else {
      f = mult(C(Math.cos(fr), Math.sin(fr)), f0)
    }

    return {
      ...z,
      ...f
    }
  }
}

function julian(pow, dist) {
  return z => {
    let a = (Math.atan2(z.im, z.re) + Math.floor(pow * Math.random()) * Math.PI * 2.0) / pow,
      r = Math.pow(dot(z, z), dist / pow * 0.5);
    return {
      ...z,
      re: Math.cos(a) * r,
      im: Math.sin(a) * r
    }
  }
}

function juliaq(pow, div) {
  let ip = div / pow,
    ip2p = (2 * Math.PI) / pow;
  return z => {
    let ang = Math.atan2(z.im, z.re) * ip + Math.floor(32767 * Math.random()) * ip2p;
    let cosa = Math.cos(ang),
      sina = Math.sin(ang),
      r = Math.pow(dot(z, z), 0.5 * ip);
    return {
      ...z,
      re: cosa * r,
      im: sina * r
    }
  }
}

function mobius(a, b, c, d) {
  return z => {
    return {
      ...z,
      ...div(add(mult(a, z), b), add(mult(c, z), d))
    }
  }
}

function murl2(c, pow) {
  return z => {
    let angle = Math.atan2(z.im, z.re) * pow;
    let cosa = Math.cos(angle),
      sina = Math.sin(angle),
      r = c * Math.pow(z.re * z.re + z.im * z.im, 0.5 * pow);
    let real = r * cosa + 1,
      imag = r * sina;
    let r1 = Math.pow(real * real + imag * imag, 1.0 / pow),
      angle1 = Math.atan2(imag, real) * 2 / pow;
    let cosa1 = Math.cos(angle1),
      sina1 = Math.sin(angle1);
    let re2 = r1 * cosa1,
      im2 = r1 * sina1,
      vp = c != -1.0 ? Math.pow(c + 1, 2.0 / pow) : 0.0;
    let r2 = vp / (r1 * r1);
    return {
      ...z,
      re: (z.re * re2 + z.im * im2) * r2,
      im: (z.im * re2 - z.re * im2) * r2
    };
  }
}

function pointSymmetry(centerX, centerY, order) {
  return z => {
    let idr = Math.floor(Math.random() * order),
      dx = z.re - centerX,
      dy = z.im - centerY,
      da = (2 * Math.PI) / order;
    let angle = idr * da;
    let cosa = Math.cos(angle),
      sina = Math.sin(angle)
    return {
      ...z,
      re: centerX + dx * cosa + dy * sina,
      im: centerY + dy * cosa - dx * sina
    }
  }
}

function rotate(theta) {
  const rad = Math.PI / 180 * theta;
  const sinTheta = -Math.sin(rad);
  const cosTheta = -Math.cos(rad);
  return z => {
    return {
      ...z,
      re: -(cosTheta * z.re + sinTheta * z.im),
      im: sinTheta * z.re - cosTheta * z.im
    }
  }
}

function scale(c) {
  if (c.hasOwnProperty('n') && c.n) {
    return function(z) {
      return {
        ...z,
        re: z.re * c.re,
        im: z.im * c.re
      }
    }
  }
  return function(z) {
    return {
      ...z,
      re: z.re * c.re,
      im: z.im * c.im
    }
  }
}

function sinusoidal() {
  return z => {
    return {
      ...z,
      re: Math.sin(z.re),
      im: Math.sin(z.im)
    }
  }
}

function smartcrop(power, radius, roundstr, roundwidth, distortion, cropmode) {
  let mode = (power > 0) == (radius > 0) ? 1 : 0,
    pow = Math.abs(power);
  let alpha = Math.PI * 2 / pow;
  let roundcoeff = roundstr / Math.sin(alpha * 0.5) / pow * 2,
    wradius = Math.abs(radius),
    radial, wpower, walpha, wcoeff;
  if (pow < 2) {
    radial = 1;
    wpower = pow * Math.PI;
    walpha = 0;
    wcoeff = 0;
  } else {
    radial = 0;
    wpower = pow;
    walpha = alpha;
    wcoeff = roundcoeff;
  }
  return function(z) {
    let ang = Math.atan2(z.im, z.re),
      rad = Math.sqrt(dot(z, z));
    let wedge, xang0, xang1, xang, coeff0, coeff1, coeff, xr, angle, wwidth, rdc;
    if (radial == 1) {
      xang0 = ang / (2 * Math.PI) + 1;
      xang = (xang0 - Math.floor(xang0)) * 2 * Math.PI;
      angle = Math.floor(Math.random() * 2) != 0 ? wpower : 0;
      if ((xang > wpower) == (mode != 1)) {
        return {
          ...z,
          re: Math.cos(angle) * rad,
          im: Math.sin(angle) * rad
        }
      } else {
        return {
          ...z
        }
      }
    } else {
      xang0 = (ang + Math.PI) / walpha;
      xang1 = xang0 - Math.floor(xang0);
      xang = xang1 < 0.5 ? xang1 : 1 - xang1;
      coeff0 = 1 / Math.cos(xang * walpha);
      wwidth = roundwidth != 1 ? Math.exp(Math.log(xang * 2) * roundwidth) * roundcoeff : xang * 2 * roundcoeff;
      coeff1 = distortion == 0 ? 1 : roundstr != 0 ? Math.abs((1 - wwidth) * coeff0 + wwidth) : coeff0;
      coeff = distortion != 1 ? Math.exp(Math.log(coeff1) * distortion) : coeff1;
      xr =  coeff * wradius;
      rdc = cropmode == -1 ? rad : xr
      f = (rad > xr) == (mode == 1) ? cropmode != 0 ? multScalar({re: Math.cos(ang), im: Math.sin(ang)}, rdc) : {re: 0, im: 0} : z;
      return {
        ...z,
        ...f
      }
    }
  }
}

function smartshape(power, roundstr, roundwidth, distortion, compensation) {
  let pow = Math.max(power, 2);
  let alpha = Math.PI * 2 / pow;
  let alphacoeff = Math.tan(alpha * 0.5) * 2,
    roundcoeff = roundstr / Math.sin(alpha * 0.5) / pow * 2,
    comp = compensation <= 0 ? 0 : 1;

  return z => {
    let dang = (Math.atan2(z.im, z.re) + Math.PI) / alpha,
      rad = Math.sqrt(dot(z, z));
    let zang1 = Math.floor(dang);
    let xang1 = dang - zang1,
      xang2, zang, sign, xang;

    if (xang1 > 0.5) {
      xang2 = 1 - xang1;
      zang = zang1 + 1;
      sign = -1;
    } else {
      xang2 = xang1;
      zang = zang1;
      sign = 1;
    }
    if (comp == 1 && distortion >= 1) {
      xang = Math.atan(xang2 * alphacoeff) / alpha;
    } else {
      xang = xang2;
    }
    let coeff0 = 1 / Math.cos(xang * alpha),
      wwidth = roundwidth != 1 ? Math.exp(Math.log(xang * 2) * roundwidth) * roundcoeff : xang * 2 * roundcoeff;
    let coeff1 = distortion == 0 ? 1 : roundstr != 0 ? Math.abs((1 - wwidth) * coeff0 + wwidth) : coeff0;
    let coeff = distortion != 1 ? Math.exp(Math.log(coeff1) * distortion) : coeff1;
    let ang = (zang + sign * xang) * alpha - Math.PI;

    return {
      ...z,
      re: Math.cos(ang) * coeff * rad,
      im: Math.sin(ang) * coeff * rad
    }
  }
}

function splits(real, imaginary) {
  return z => {
    const xoff = z.re > 0 ? real : -real,
      yoff = z.im > 0 ? imaginary : -imaginary;
    return {
      ...z,
      re: z.re + xoff,
      im: z.im + yoff
    }
  }
}

function tileHelp() {
  return z => {
    let x = z.re;
    let val = Math.cos((x > 0 ? x - Math.floor(x) : x + Math.floor(-x)) * Math.PI),
      fpx;
    if (val < Math.random() * 2 - 1) {
      fpx = x > 0 ? -1 : 1
    } else {
      fpx = 0
    }
    return {
      ...z,
      re: z.re + fpx,
      im: z.im
    }
  }
}

function tileLog(spread) {
  return z => {
    return {
      ...z,
      re: z.re + Math.floor(Math.log(Math.random()) * (Math.random() < 0.5 ? spread : -spread) + 0.5),
      im: z.im
    }
  }
}

function translate(real, imaginary) {
  return z => {
    return {
      ...z,
      re: z.re + real,
      im: z.im + imaginary
    }
  }
}

function trigCosh() {
  return z => {
    return {
      ...z,
      ...cosh(z)
    }
  }
}

function trigExp() {
  return z => {
    return {
      ...z,
      ...exp(z)
    }
  }
}

function trigLog() {
  return z => {
    return {
      ...z,
      ...log(z)
    }
  }
}

function trigSinh() {
  return z => {
    return {
      ...z,
      ...sinh(z)
    }
  }
}

function trigTanh() {
  return z => {
    return {
      ...z,
      ...tanh(z)
    }
  }
}

function unbubble() {
  return z => {
    let r = dot(z, z);
    let b = (Math.SQRT2 - Math.sqrt(2 - r)) / r;
    return {
      ...z,
      re: z.re * b,
      im: z.im * b
    }
  }
}

/* colors */

function color(color) {
  let col = color;
  if (color.hasOwnProperty('h')) {
    col = hslToRgb(color.h, color.s, color.l);
  } else if (!color.hasOwnProperty('red')) {
    col = colorRGB(...arguments);
  }
  return z => {
    return {
      ...z,
      red: col.red,
      green: col.green,
      blue: col.blue
    }
  }
}

function lerp(colorA, colorB, weight = 0.5) {
  weight = Math.max(0, Math.min(1, weight));
  if (colorA.hasOwnProperty('h')) {
    if (colorB.hasOwnProperty('h')) {
      return hslToRgb(
        Math.abs((colorA.h%1) - (colorB.h%1)) > 0.5 ?
          (((colorA.h+(colorA.h%1>colorB.h?0:1)) * (1 - weight) + (colorB.h+(colorA.h%1>colorB.h?1:0)) * weight + 1)+1) % 1 :
          colorA.h * (1 - weight) + colorB.h * weight,
        colorA.s * (1 - weight) + colorB.s * weight,
        colorA.l * (1 - weight) + colorB.l * weight);
    } else {
      colorB = rgbToHsl(colorB.red, colorB.green, colorB.blue);
      return hslToRgb(
        Math.abs((colorA.h%1) - (colorB.h%1)) > 0.5 ?
          (((colorA.h+(colorA.h%1>colorB.h?0:1)) * (1 - weight) + (colorB.h+(colorA.h%1>colorB.h?1:0)) * weight + 1)+1) % 1 :
          colorA.h * (1 - weight) + colorB.h * weight,
        colorA.s * (1 - weight) + colorB.s * weight,
        colorA.l * (1 - weight) + colorB.l * weight);
    }
  } else {
    if (colorB.hasOwnProperty('h')) {
      colorA = rgbToHsl(colorA.red, colorA.green, colorA.blue);
      return hslToRgb(
        Math.abs((colorA.h%1) - (colorB.h%1)) > 0.5 ?
          (((colorA.h+(colorA.h%1>colorB.h?0:1)) * (1 - weight) + (colorB.h+(colorA.h%1>colorB.h?1:0)) * weight + 1)+1) % 1 :
          colorA.h * (1 - weight) + colorB.h * weight,
        colorA.s * (1 - weight) + colorB.s * weight,
        colorA.l * (1 - weight) + colorB.l * weight);
    } else {
      return lerpRGB(colorB.red, colorB.green, colorB.blue, weight)(colorA);
    }
  }
}

function lerpColor(color, weight) {
  return z => {
    return {
      ...z,
      ...lerp(z, color, weight)
    }
  }
}

function colorRGB(r, g, b) {
  return {
    red: r,
    green: g,
    blue: b
  }
}

function colorHSL(h, s, l) {
  return {
    h: h,
    s: s,
    l: l
  };
}

function hue2rgb(p, q, t) {
  if (t < 0) t += 1;
  if (t > 1) t -= 1;
  if (t < 1 / 6) return p + (q - p) * 6 * t;
  if (t < 1 / 2) return q;
  if (t < 2 / 3) return p + (q - p) * (2 / 3 - t) * 6;
  return p;
}

function hslToRgb(h, s, l) {
  let r, g, b;

  if (s == 0) {
    r = g = b = l; // achromatic
  } else {
    let q = l < 0.5 ? l * (1 + s) : l + s - l * s;
    let p = 2 * l - q;

    r = hue2rgb(p, q, h + 1 / 3);
    g = hue2rgb(p, q, h);
    b = hue2rgb(p, q, h - 1 / 3);
  }

  return {
    red: r,
    green: g,
    blue: b
  };
}

function rgbToHsl(r, g, b) {
  let max = Math.max(r, g, b),
    min = Math.min(r, g, b);
  let h, s, l = (max + min) / 2;

  if (max == min) {
    h = s = 0; // achromatic
  } else {
    let d = max - min;
    s = l > 0.5 ? d / (2 - max - min) : d / (max + min);

    switch (max) {
      case r:
        h = (g - b) / d + (g < b ? 6 : 0);
        break;
      case g:
        h = (b - r) / d + 2;
        break;
      case b:
        h = (r - g) / d + 4;
        break;
    }

    h /= 6;
  }

  return {
    h: h,
    s: s,
    l: l
  };
}

function hslShift(h, s, l) {
  return z => {
    let hsl = rgbToHsl(z.red, z.green, z.blue);
    return {
      ...z,
      ...hslToRgb((hsl.h + h + 1) % 1, hsl.s + s, hsl.l + l)
    }
  }
}

function setHue(hue) {
  return z => {
    let hsl = rgbToHsl(z.red, z.green, z.blue);
    return {
      ...z,
      ...hslToRgb(hue, hsl.s, hsl.l)
    }
  }
}

function setSaturation(saturation) {
  return z => {
    let hsl = rgbToHsl(z.red, z.green, z.blue);
    return {
      ...z,
      ...hslToRgb(hsl.h, saturation, hsl.l)
    }
  }
}

function setLightness(lightness) {
  return z => {
    let hsl = rgbToHsl(z.red, z.green, z.blue);
    return {
      ...z,
      ...hslToRgb(hsl.h, hsl.s, lightness)
    }
  }
}

function lerpRGB(r, g, b, weight = 0.5) {
  weight = Math.max(0, Math.min(1, weight));
  return z => {
    return {
      ...z,
      ...{
        red: (z.red * (1 - weight) + r * weight),
        green: (z.green * (1 - weight) + g * weight),
        blue: (z.blue * (1 - weight) + b * weight)
      }
    }
  }
}

function lerpHSL(h, s, l, weight = 0.5) {
  weight = Math.max(0, Math.min(1, weight));
  return z => {
    let hsl = rgbToHsl(z.red, z.green, z.blue);
    return {
      ...z,
      ...hslToRgb(
        Math.abs((hsl.h%1) - (h%1)) > 0.5 ?
          (((hsl.h+(hsl.h%1>h%1?0:1)) * (1 - weight) + (h+(hsl.h%1>h%1?1:0)) * weight)+1)%1:
          (hsl.h%1) * (1 - weight) + (h%1) * weight,
        hsl.s * (1 - weight) + s * weight,
        hsl.l * (1 - weight) + l * weight
      )
    }
  }
}

function normalizeRGB(color) {
  let magnitude = Math.max(color.red, color.green, color.blue);
  return {
    red: color.red / magnitude,
    green: color.green / magnitude,
    blue: color.blue / magnitude
  }
}

function normalizeColor() {
  return z => {
    let magnitude = Math.max(z.red, z.green, z.blue);
    return {
      ...z,
      red: z.red / magnitude,
      green: z.green / magnitude,
      blue: z.blue / magnitude
    }
  }
}

function brightenRGB(color, amount) {
  return {
    red: Math.min(1, color.red * amount),
    green:  Math.min(1, color.green * amount),
    blue:  Math.min(1, color.blue * amount)
  }
}

function brighten(amount) {
  return z => {
    return {
      ...z,
      alpha: z.alpha * amount,
    }
  }
}

function repeatingGradient(colors) {
  if (colors.length < 1) {
    throw "not enough colors";
  }
  if (colors.length === 1) {
    let col = colors[0];
    if (colors[0].hasOwnProperty('h')) {
      col = hslToRgb(colors[0].h, colors[0].s, colors[0].l);
    }
    return z => {
      return {
        ...z,
        ...col
      }
    }
  }
  return z => {
    let at = (z.im - Math.floor(z.im)) * colors.length;
    return {
      ...z,
      ...lerp(colors[Math.floor(at) % colors.length], colors[Math.ceil(at) % colors.length], at % 1)
    }
  }
}


function gradient(colors) {
  if (colors.length < 1) {
    throw "not enough colors";
  }
  if (colors.length === 1) {
    let col = colors[0];
    if (colors[0].hasOwnProperty('h')) {
      col = hslToRgb(colors[0].h, colors[0].s, colors[0].l);
    }
    return z => {
      return {
        ...z,
        ...col
      }
    }
  }
  return z => {
    let at = Math.max(0, Math.min(1, z.im + 0.5)) * (colors.length - 1);
    return {
      ...z,
      ...lerp(colors[Math.floor(at) % colors.length], colors[Math.ceil(at) % colors.length], at % 1)
    }
  }
}

function gamma(gamma) {
  return z => {
    return {
      ...z,
      red: Math.pow(z.red, 1 / gamma),
      green: Math.pow(z.green, 1 / gamma),
      blue: Math.pow(z.blue, 1 / gamma)
    }
  }
}

/* description s*/

const BUILT_IN_TRANSFORMS = {
  reset: reset,
  arcsinh: arcsinh,
  arctanh: arctanh,
  blurCircle: blurCircle,
  blurGasket: blurGasket,
  blurGaussian: blurGaussian,
  blurSine: blurSine,
  blurSquare: blurSquare,
  bubble: bubble,
  circleInv: circleInv,
  cpow: cpow,
  cylinder: cylinder,
  disc: disc,
  hypershape: hypershape,
  hypershift: hypershift,
  hypertile3: hypertile3,
  julian: julian,
  juliaq: juliaq,
  mobius: mobius,
  murl2: murl2,
  pointSymmetry: pointSymmetry,
  rotate: rotate,
  scale: scale,
  sinusoidal: sinusoidal,
  smartcrop: smartcrop,
  smartshape: smartshape,
  splits: splits,
  tileHelp: tileHelp,
  tileLog: tileLog,
  translate: translate,
  trigCosh: trigCosh,
  trigExp: trigExp,
  trigLog: trigLog,
  trigSinh: trigSinh,
  trigTanh: trigTanh,
  unbubble: unbubble,
  //color transfomrs
  brighten: brighten,
  color: color,
  gamma: gamma,
  gradient: gradient,
  repeatingGradient: repeatingGradient,
  hslShift: hslShift,
  lerpColor: lerpColor,
  lerpHSL: lerpHSL,
  lerpRGB: lerpRGB,
  normalizeColor: normalizeColor,
  setHue: setHue,
  setSaturation: setSaturation,
  setLightness: setLightness,
};

const BUILT_IN_TRANSFORMS_PARAMS = {
  reset: [],
  arcsinh: [],
  arctanh: [],
  blurCircle: [],
  blurGasket: [],
  blurGaussian: [{
    name: "pow",
    type: "number",
    default: 1
  }],
  blurSine: [{
    name: "pow",
    type: "number",
    default: 1
  }],
  blurSquare: [],
  bubble: [],
  circleInv: [],
  cpow: [{
      name: "rotation",
      type: "number",
      default: 1
    },
    {
      name: "depth",
      type: "number",
      default: 1
    },
    {
      name: "divisor",
      type: "number",
      default: 1
    },
    {
      name: "spread",
      type: "number",
      default: 1
    }
  ],
  cylinder: [],
  disc: [],
  hypershape: [{
    name: "pow",
    type: "number",
    default: 3
  }],
  hypershift: [{
    name: "p",
    type: "complex",
    default: {
      re: 0,
      im: 0
    }
  }],
  hypertile3: [{
      name: "p",
      type: "number",
      default: 3
    },
    {
      name: "q",
      type: "number",
      default: 3
    },
    {
      name: "r",
      type: "number",
      default: 4
    },
    {
      name: "d",
      type: "number",
      default: 0.5
    }
  ],
  julian: [{
      name: "pow",
      type: "number",
      default: 1
    },
    {
      name: "dist",
      type: "number",
      default: 1
    }
  ],
  juliaq: [{
      name: "pow",
      type: "number",
      default: 1
    },
    {
      name: "div",
      type: "number",
      default: 1
    }
  ],
  mobius: [{
      name: "a",
      type: "complex",
      default: {
        re: 1,
        im: 0
      }
    },
    {
      name: "b",
      type: "complex",
      default: {
        re: 0,
        im: 0
      }
    },
    {
      name: "c",
      type: "complex",
      default: {
        re: 0,
        im: 0
      }
    },
    {
      name: "d",
      type: "complex",
      default: {
        re: 1,
        im: 0
      }
    }
  ],
  murl2: [{
      name: "c",
      type: "number",
      default: 0
    },
    {
      name: "pow",
      type: "number",
      default: 2
    }
  ],
  pointSymmetry: [{
      name: "centerX",
      type: "number",
      default: 0
    },
    {
      name: "centerY",
      type: "number",
      default: 0
    },
    {
      name: "order",
      type: "number",
      default: 1
    }
  ],
  rotate: [{
    name: "theta",
    type: "number",
    default: 0
  }],
  scale: [{
    name: "c",
    type: "complex",
    default: {
      re: 1,
      im: 1
    }
  }],
  sinusoidal: [],
  smartcrop: [{
      name: "power",
      type: "number",
      default: 4
    },
    {
      name: "radius",
      type: "number",
      default: 1
    },
    {
      name: "roundstr",
      type: "number",
      default: 0
    },
    {
      name: "roundwidth",
      type: "number",
      default: 1
    },
    {
      name: "distortion",
      type: "number",
      default: 1
    },
    {
      name: "cropmode",
      type: "number",
      default: 2
    }
  ],
  smartshape: [{
      name: "power",
      type: "number",
      default: 4
    },
    {
      name: "roundstr",
      type: "number",
      default: 0
    },
    {
      name: "roundwidth",
      type: "number",
      default: 1
    },
    {
      name: "distortion",
      type: "number",
      default: 1
    },
    {
      name: "compensation",
      type: "number",
      default: 1
    }
  ],
  splits: [{
      name: "real",
      type: "number",
      default: 0
    },
    {
      name: "imaginary",
      type: "number",
      default: 0
    }
  ],
  tileHelp: [],
  tileLog: [{
    name: "spread",
    type: "number",
    default: 1
  }],
  translate: [{
      name: "real",
      type: "number",
      default: 0
    },
    {
      name: "imaginary",
      type: "number",
      default: 0
    }
  ],
  trigCosh: [],
  trigExp: [],
  trigLog: [],
  trigSinh: [],
  trigTanh: [],
  unbubble: [],
  //color transforms
  brighten: [{
    name: "amount",
    type: "number",
    default: 1
  }],
  color: [{
    name: "color",
    type: "object",
    default: {
      h: 0,
      s: 1,
      l: 0.5
    }
  }],
  gamma: [{
    name: "gamma",
    type: "number",
    default: "2.2"
  }],
  gradient: [{
    name: "colorA",
    type: "array",
    default: [{
      red: 1,
      green: 1,
      blue: 1
    }, {
      red: 1,
      green: 0,
      blue: 0
    }]
  }],
  repeatingGradient: [{
    name: "colorA",
    type: "array",
    default: [{
      red: 1,
      green: 1,
      blue: 1
    }, {
      red: 1,
      green: 0,
      blue: 0
    }]
  }],
  hslShift: [{
    name: "h",
    type: "number",
    default: 0
  }, {
    name: "s",
    type: "number",
    default: 0
  }, {
    name: "l",
    type: "number",
    default: 0
  }],
  lerpHSL: [{
      name: "h",
      type: "number",
      default: 0
    },
    {
      name: "s",
      type: "number",
      default: 1
    },
    {
      name: "l",
      type: "number",
      default: 0.5
    }, {
      name: "weight",
      type: "number",
      default: 0.5
    }
  ],
  lerpColor: [{
    name: "color",
    type: "object",
    default: {
      red: 1,
      green: 1,
      blue: 1
    }
  }, {
    name: "weight",
    type: "number",
    default: 0.5
  }],
  lerpRGB: [{
      name: "r",
      type: "number",
      default: 1
    },
    {
      name: "g",
      type: "number",
      default: 1
    },
    {
      name: "b",
      type: "number",
      default: 1
    }, {
      name: "weight",
      type: "number",
      default: 0.5
    }
  ],
  normalizeColor: [],
  setHue: [{
    name: "hue",
    type: "number",
    default: 0
  }],
  setSaturation: [{
    name: "saturation",
    type: "number",
    default: 1
  }],
  setLightness: [{
    name: "lightness",
    type: "number",
    default: 0.5
  }],
};
