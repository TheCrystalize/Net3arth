/*helper functions*/
function C(real, imaginary) {
  return {
    re: real,
    im: imaginary
  };
}

function div(z, c) {
  const s = 1 / (c.re * c.re + c.im * c.im);
  return {
    re: (z.re * c.re + z.im * c.im) * s,
    im: (z.im * c.re - z.re * c.im) * s
  }
}

function add(z, c) {
  return {
    re: z.re + c.re,
    im: z.im + c.im
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
  const s = Math.sqrt(z.re * z.re + z.im * z.im),
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

/*Transforms*/
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

function mobius(a, b, c, d) {
  return function(z) {
    return {
      ...z,
      ...div(add(mult(a, z), b), add(mult(c, z), d))
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

function translate(real, imaginary) {
  return function(z) {
    return {
      ...z,
      re: z.re + real,
      im: z.im + imaginary
    }
  }
}

function rotate(theta) {
  const sinTheta = -Math.sin(theta);
  const cosTheta = -Math.cos(theta);
  return function(z) {
    return {
      ...z,
      re: cosTheta * z.re + sinTheta * z.im,
      im: sinTheta * z.re - cosTheta * z.im
    }
  }
}

function blurCircle() {
  return function(z) {
    let a = Math.random() * Math.PI * 2,
      r = Math.sqrt(Math.random());
    let ans = {
      re: Math.cos(a) * r,
      im: Math.sin(a) * r
    }
    return {
      ...z,
      ...ans
    }
  }
}

function hypertile3(p, q, r, shift) {
  let rad = Math.pi / 180,
    o = Math.acosh((Math.cos(Math.PI / p) + Math.cos(Math.PI / q) * Math.cos(Math.PI / r)) / (Math.sin(Math.PI / q) * Math.sin(Math.PI / r))),
    a = Math.asinh(Math.sin(Math.PI / q) / Math.sin(Math.PI / p) * Math.sinh(o)),
    b = Math.asinh(Math.sin(Math.PI / r) * Math.sinh(a)),
    c = Math.asinh(Math.sin(Math.PI / r) * Math.sinh(a) / Math.sin(Math.PI / q)),
    h = Math.tanh(b / 2),
    b1 = Math.tanh(Math.acosh(Math.cosh(c) / Math.cosh(b)) / 2),
    b2 = Math.tanh(Math.acosh(Math.cosh(a) / Math.cosh(b)) / 2),
    rot1 = 360 / p * rad,
    rot2 = 360 / q * rad,
    rot3 = 360 / r * rad,
    r01 = C(Math.cos(rot1), Math.sin(rot1)),
    r02 = C(Math.cos(rot2), Math.sin(rot2)),
    r03 = C(Math.cos(rot3), Math.sin(rot3));

  const c1 = C(h, 0),
    c2 = C(0, b1),
    c3 = C(0, b2);

  return function(z) {
    let n, prf, z0;
    if(shift < 0.25) {
      n = 0
      pfr = 0
      z0 = z
    } else if(shift < 0.5) {
      n = p
      pfr = rot1
      z0 = div(add(z, c1), add(mult(c1, z), 1))
    } else if(shift < 0.75) {
      n = q
      pfr = rot2
      z0 = div(add(z, -c2), add(mult(c2, z), 1))
    } else {
      n = r
      pfr = rot3
      z0 = div(add(z, c3), add(mult(-c3, z), 1))
    }

    let m0 = div(add(z0, -c1), add(mult(-c1, z0), 1)),
      r1 = r01 * m0,
      m0f = div(add(r1, c1), add(mult(c1, r1), 1)),
      m1 = div(add(z0, c2), add(mult(-c2, z0), 1)),
      r2 = r02 * m1,
      m1f = div(add(r2, -c2), add(mult(c2, r2), 1)),

      m2 = div(add(z0, -c3), add(mult(c3, z0), 1)),
      r3 = r03 * m2,
      m2f = div(add(r3, c3), add(mult(-c3, r3), 1)),

      fr = Math.floor(Math.random() * n) * pfr,
      rnd = Math.random(),
      f3, f0, f;
    if(rnd < 1 / 3) {
      f3 = m0f
    } else if(rnd < 2 / 3) {
      f3 = m1f
    } else {
      f3 = m2f
    };
    if(shift < 0.25) {
      f0 = f3
    } else if(shift < 0.5) {
      f0 = div(add(f3, -c1), add(mult(-c1, f3), 1))
    } else if(shift < 0.75) {
      f0 = div(add(f3, c2), add(mult(-c2, f3), 1))
    } else {
      f0 = div(add(f3, -c3), add(mult(c3, f3), 1))
    };
    if(shift < 0.25) {
      f = f0
    } else {
      f = C(Math.cos(fr), Math.sin(fr)) * f0
    };
    let ans = f;
    return {
      ...z,
      ...ans
    }
  }
}

const BUILT_IN_TRANSFORMS = {
  arcsinh: arcsinh,
  splits: splits,
  mobius: mobius,
  scale: scale,
  translate: translate,
  rotate: rotate,
  blurCircle: blurCircle,
  hypertile3: hypertile3
};

const BUILT_IN_TRANSFORMS_PARAMS = {
    arcsinh: [],
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
    scale: [{
      name: "s",
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
    rotate: [{
      name: "theta",
      type: "number",
      default: 0
    }],
    blurCircle: [],
    hypertile3: [
      {
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
    ]
  };
