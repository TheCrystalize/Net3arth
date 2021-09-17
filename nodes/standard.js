function arcsinh() {
  return function(z) {
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
  return function(z) {
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
  return function(z) {
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
  return function(z) {
    let a = Math.random() * Math.PI * 2,
      r = 1 / Math.sqrt(Math.random() - Math.random());
    return {
      ...z,
      re: Math.random() - 0.5,
      im: Math.sin(a) * r
    }
  }
}

function blurSine(pow) {
  return function(z) {
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
  return function(z) {
    return {
      ...z,
      re: Math.random() - 0.5,
      im: Math.random() - 0.5
    }
  }
}

function circleInv() {
  return function(z) {
    let r = 1 / dot(z, z)
    return {
      ...z,
      re: z.re * r,
      im: -z.im * r
    }
  }
}

function hypershift(p) {
  return function(z) {
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

  return function(z) {
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
  return function(z) {
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
  return function(z) {
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
  return function(z) {
    return {
      ...z,
      ...div(add(mult(a, z), b), add(mult(c, z), d))
    }
  }
}

function murl2(c, pow) {
  return function(z) {
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
  return function(z) {
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
  const th = theta * Math.PI / 2;
  const sinTheta = -Math.sin(th);
  const cosTheta = -Math.cos(th);
  return function(z) {
    return {
      ...z,
      re: cosTheta * z.re + sinTheta * z.im,
      im: sinTheta * z.re - cosTheta * z.im
    }
  }
}

function scale(s) {
  return function(z) {
    return {
      ...z,
      ...multScalar(z, s)
    }
  }
}

function splits(x, y) {
  return function(z) {
    const xoff = z.re > 0 ? x : -x,
      yoff = z.im > 0 ? y : -y;
    return {
      ...z,
      re: z.re + xoff,
      im: z.im + yoff
    }
  }
}

function tileHelp(width) {
  return function(z) {
    let x = z.re / width;
    let val = Math.cos((x > 0 ? x - Math.floor(x) : x + Math.floor(x)) * Math.PI),
      fpx;
    if (val < Math.random() * 2 - 1) {
      fpx = x > 0 ? -width : width
    } else {
      fpx = 0
    }
    return {
      ...z,
      re: z.re + fpx
    }
  }
}

function tileLog(spread) {
  return function(z) {
    return {
      ...z,
      re: z.re + Math.floor(Math.log(Math.random()) * (Math.random() < 0.5 ? spread : -spread) + 0.5)
    }
  }
}

function translate(real, imaginary) {
  return function(z) {
    return {
      ...z,
      re: z.re + real,
      im: z.im + imaginary
    }
  }
}

function trigCosh() {
  return function(z) {
    return {
      ...z,
      ...cosh(z)
    }
  }
}

function trigExp() {
  return function(z) {
    return {
      ...z,
      ...exp(z)
    }
  }
}

function trigLog() {
  return function(z) {
    return {
      ...z,
      ...log(z)
    }
  }
}

function trigSinh() {
  return function(z) {
    return {
      ...z,
      ...sinh(z)
    }
  }
}

function trigTanh() {
  return function(z) {
    return {
      ...z,
      ...tanh(z)
    }
  }
}

function unbubble() {
  return function(z) {
    let r = dot(z, z);
    let b = (Math.SQRT2 - Math.sqrt(2 - r)) / r;
    return {
      ...z,
      re: z.re * b,
      im: z.im * b
    }
  }
}

const BUILT_IN_TRANSFORMS = {
  arcsinh: arcsinh,
  arctanh: arctanh,
  blurCircle: blurCircle,
  blurGasket: blurGasket,
  blurSine: blurSine,
  blurSquare: blurSquare,
  circleInv: circleInv,
  hypershift: hypershift,
  hypertile3: hypertile3,
  julian: julian,
  juliaq: juliaq,
  mobius: mobius,
  murl2: murl2,
  pointSymmetry: pointSymmetry,
  rotate: rotate,
  scale: scale,
  splits: splits,
  tileHelp: tileHelp,
  tileLog: tileLog,
  translate: translate,
  trigCosh: trigCosh,
  trigExp: trigExp,
  trigLog: trigLog,
  trigSinh: trigSinh,
  trigTanh: trigTanh,
  unbubble: unbubble
};

const BUILT_IN_TRANSFORMS_PARAMS = {
  arcsinh: [],
  arctanh: [],
  blurCircle: [],
  blurGasket: [],
  blurSine: [{
    name: "pow",
    type: "number",
    default: 1
  }],
  blurSquare: [],
  circleInv: [],
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
    name: "s",
    type: "number",
    default: 1
  }],
  splits: [{
      name: "x",
      type: "number",
      default: 0
    },
    {
      name: "y",
      type: "number",
      default: 0
    }
  ],
  tileHelp: [{
    name: "width",
    type: "number",
    default: 1
  }],
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
  unbubble: []
};
