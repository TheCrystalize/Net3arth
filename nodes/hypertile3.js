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

const hypertile3Data = {
  name: "hypertile3",
  parameters: [
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
}
