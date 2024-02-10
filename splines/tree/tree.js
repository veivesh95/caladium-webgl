function rami(b1, b2, th1, th2, nfi, R0, R1, R2, l0, l1, l2, a, stu) {
  var bc = new Float32Array(12), // cpc, bnc, bnc2, bn0top
    b0 = new Float32Array(b1),
    ssu = (su - 1) / stu + 1;
  var co = Math.cos(th1),
    si = Math.sin(th1),
    ta2 = Math.tan(th1 * 0.5),
    fi = (Math.PI * nfi) / (ssu - 1),
    cof = Math.cos(fi),
    sif = Math.sin(fi);
  for (var i = 0; i < 3; i++) {
    b0[i + 6] = sif * b1[i + 3] + cof * b1[i + 6]; // fi rot
    b0[i + 3] = b1[i + 3] = b2[i + 3] = cof * b1[i + 3] - sif * b1[i + 6];
    bc[i] = b1[i + 9] + l0 * b1[i];
    bc[i + 3] = ta2 * b0[i] + b0[i + 6]; // to look like an even stem
    b1[i + 6] = si * b0[i] + co * b0[i + 6]; // th rot
    b1[i] = co * b0[i] - si * b0[i + 6];
    b1[i + 9] = bc[i] + l1 * b1[i];
  }
  b1[12] -= nfi / (ssu - 1);
  b2[12] = b1[12];
  b2[13] = b1[13];
  var Rm = (R0 * l0 + R1 * l1) / (l0 + l1);
  patch(b1, R0, R1, Rm, 0, 0, l0 + l1); // side 1
  var th = 0.5 * (th1 + th2),
    dth = 0.5 * (th1 - th2),
    si2 = Math.sin(dth * (1 - a)),
    t = th + dth * a;
  co = Math.cos(t);
  si = Math.sin(t);
  for (var i = 0; i < 3; i++) bc[i + 3] = (b0[i] * co - b0[i + 6] * si) / si2; // -bnc1
  t = th - dth * a;
  co = Math.cos(t);
  si = Math.sin(t);
  var Rm2 = (R0 * l0 + R2 * l2) / (l0 + l2);
  for (var i = 0; i < 3; i++) {
    bc[i + 6] = (b0[i] * co - b0[i + 6] * si) / -si2; // bnc2
    bc[i + 9] = -b0[i] * l0 - (-bc[i + 3] * Rm + bc[i + 6] * Rm2) * 0.5;
  }
  patch(b1, R0, R1, Rm, su - 1, 1, l0 + l1); // top 1
  co = Math.cos(th2);
  si = Math.sin(th2);
  ta2 = Math.tan(th2 * 0.5);
  for (var i = 0; i < 3; i++) {
    bc[i + 3] = ta2 * b0[i] + b0[i + 6];
    b2[i + 6] = si * b0[i] + co * b0[i + 6]; // th rot
    b2[i] = co * b0[i] - si * b0[i + 6];
    b2[i + 9] = bc[i] + l2 * b2[i];
  }
  patch(b2, R0, R2, Rm2, su - 1, 0, l0 + l2); // side 2
  for (var i = 0; i < 3; i++) {
    bc[i + 3] = bc[i + 6];
    bc[i + 9] = -bc[i + 9];
  }
  patch(b2, R0, R2, Rm2, 0, 1, l0 + l2); // top 2
  b1[13] += 0.5 * (l0 + l1);
  b2[13] += 0.5 * (l0 + l2);

  function patch(b, R0, R1, Rm, u0, top, ls) {
    ls *= 0.5;
    var t = off / 8,
      tm = t + ssu,
      to,
      bo,
      si,
      co,
      c0,
      c1;
    for (var j = 0; j < sv - 1; j++) {
      for (var i = 0; i < su - 1; i += stu) {
        ind[pi++] = t++;
        ind[pi++] = t;
        ind[pi++] = tm;
        ind[pi++] = tm++;
        ind[pi++] = t;
        ind[pi++] = tm;
      }
      t++;
      tm++;
    }
    if (top) {
      c0 = u0 ? -b0[0] : b0[0];
      c1 = (c0 + b[6]) * 0.5;
    } else {
      c0 = b0[6];
      c1 = bc[3];
    }
    for (var v = 0; v < sv; v++) {
      var t = v / (sv - 1),
        t1 = 1 - t,
        B0 = t1 * t1,
        B1 = 2 * t * t1,
        B2 = t * t;
      for (var u = 0; u < su; u += stu) {
        si = sin[u + u0];
        co = cos[u + u0];
        for (var i = 0; i < 3; i++) {
          bo = b[i + 3] * co;
          to = top ? bc[i + 9] * si : b0[i + 6] * si * R0;
          pt[off++] =
            (b0[i + 9] + to + bo * R0) * B0 +
            (bc[i] + (bc[i + 3] * si + bo) * Rm) * B1 +
            (b[i + 9] + (b[i + 6] * si + bo) * R1) * B2;
        }
        var li = 0.5 * (b[3] * co + (B0 * c0 + B1 * c1 + B2 * b[6]) * si);
        if (li < 0) li = 0;
        pt[off++] = 0.6 + li;
        pt[off++] = 0.3 + li;
        pt[off++] = 0;
        pt[off++] = ls * t + b[13];
        pt[off++] = u / (su - 1) + b[12];
      }
    }
  } // end patch
}
function base(b, nfi, R0, R2, dR, l0, l1, redu) {
  var t = off / 8,
    tm = t + sv;
  for (var u = 0; u < 2 * su - 2; u++) {
    for (var j = 0; j < sv - 1; j++) {
      ind[pi++] = t++;
      ind[pi++] = tm;
      ind[pi++] = t;
      ind[pi++] = tm++;
      ind[pi++] = tm;
      ind[pi++] = t;
    }
    t++;
    tm++;
  }
  var fi = (Math.PI * nfi) / (su - 1),
    cof = Math.cos(fi),
    sif = Math.sin(fi);
  var cp = new Float32Array(6);
  for (var i = 0; i < 3; i++) {
    cp[i] = b[i + 9];
    cp[i + 3] = b[i + 9] + l0 * b[i];
    t = sif * b[i + 3] + cof * b[i + 6]; // fi rot
    b[i + 3] = cof * b[i + 3] - sif * b[i + 6];
    b[i + 6] = t;
    b[i + 9] = cp[i + 3] + l1 * b[i];
  }
  b[12] -= nfi / (su - 1);
  var du = (14 * Math.PI) / (2 * su - 2);
  for (var u = 0; u < 2 * su - 1; u++) {
    var si = sin[u],
      co = cos[u],
      si2 = si,
      co2 = co;
    if (redu && u & 1) {
      si2 = (sin[u - 1] + sin[u + 1]) * 0.5;
      co2 = (cos[u - 1] + cos[u + 1]) * 0.5;
    }
    var Ru = R0 + dR * Math.sin(du * u);
    for (var v = 0; v < sv; v++) {
      var t = v / (sv - 1),
        t1 = 1 - t,
        B0 = t1 * t1,
        B1 = 2 * t * t1,
        B2 = t * t;
      for (var i = 0; i < 3; i++) {
        pt[off++] =
          (cp[i] + (b[i + 6] * si + b[i + 3] * co) * Ru) * B0 +
          (cp[i + 3] + (b[i + 6] * si + b[i + 3] * co) * R2) * B1 +
          (b[i + 9] + (b[i + 6] * si2 + b[i + 3] * co2) * R2) * B2;
      }
      var li = 0.5 * (b[3] * co + b[6] * si);
      if (li < 0) li = 0;
      pt[off++] = 0.6 + li;
      pt[off++] = 0. + li;
      pt[off++] = 0;
      pt[off++] = t + b[13];
      pt[off++] = u / (su - 1) + b[12];
    }
  }
  b[13]++;
}

function leaves3(b, R, green, age) {
  if (b[1] > 0) {
    b[1] = b[1] * (1 - b[1] * b[1]);
    var no = 0;
    for (var j = 0; j < 3; j++) no += b[j] * b[j];
    no = 1 / Math.sqrt(no);
    for (var j = 0; j < 3; j++) b[j] *= no;
  }

  leaf3(b, R, 1, 0.1, 0, green, age); // change stem to leaf position
  // console.log('green', green)
  console.log("R", R);
  console.log("b", b);

  // leaf3(b, R, 0.7, -0.1, 0.9, green);
  // leaf3(b, R, 0.7, 0.1, -1.8, green);
  // leaf3(b, R, 1.4, -0.15, 0.4, green);
  // leaf3(b, R, 1.5, -0.05, 0.9, green);
}
function leaf30(b, R, lb, y0, ang, green, age) {
  // autumn
  console.log(age);
  var dis = Math.round(b[13] + 2);
  var t0 = off / 8,
    r1 = 0.6 * age,
    rm = 0.5 * age + 1,
    r2 = 0.6 * age,
    si = Math.sin(ang),
    co = Math.cos(ang),
    fi;
  var t = b[0] * co - b[2] * si;
  b[2] = b[0] * si + b[2] * co;
  b[0] = t;
  b[1] += y0;
  var bn = new Float32Array(6);
  var no = Math.sqrt(b[0] * b[0] + b[2] * b[2]);
  bn[0] = b[2] / no;
  bn[2] = -b[0] / no;
  bn[4] = 0.2 + age * 0.15; // to change the depth of the leaf
  for (var j = 0; j < 3; j++) pt[off++] = b[j + 9] + lb * b[j];

  // for the color
  var li = 0.7 + 0.2 * (1 - Math.abs(b[1])),
    rli = 0;
  if (!green) {
    rli = 0.4 + li * rnd[irnd++];
    li = 0.3 + li * rnd[irnd] * rnd[irnd++];
  }

  pt[off++] = 0.7 * rli;
  pt[off++] = 0.6 * li;
  pt[off++] = 0;
  pt[off++] = dis;
  pt[off++] = 0.5;
  t = t0 + 1;
  var sd = Math.sin(Math.PI * 0.1),
    cd = Math.cos(Math.PI * 0.1),
    si = 0,
    co = -1,
    tmp;
  for (var i = 0; i < 5; i++) {
    for (var j = 0; j < 3; j++)
      pt[off++] =
        b[j + 9] +
        lb * b[j] +
        bn[j + 3] +
        b[j] * (r1 * co + 0.1) +
        bn[j] * r1 * si;
    pt[off++] = rli;
    pt[off++] = li;
    pt[off++] = 0;
    pt[off++] = 0.4 + dis;
    pt[off++] = 1;
    tmp = si * cd + co * sd;
    co = co * cd - si * sd;
    si = tmp;
    for (var j = 0; j < 3; j++)
      pt[off++] =
        b[j + 9] +
        lb * b[j] +
        bn[j + 3] +
        b[j] * (rm * co + 0.15) +
        bn[j] * rm * si;
    pt[off++] = rli;
    pt[off++] = li;
    pt[off++] = 0;
    pt[off++] = 0.75 + dis;
    pt[off++] = 1;
    tmp = si * cd + co * sd;
    co = co * cd - si * sd;
    si = tmp;
    for (var j = 0; j < 3; j++)
      pt[off++] =
        b[j + 9] +
        lb * b[j] +
        bn[j + 3] +
        b[j] * (r2 * co + 0.2) +
        bn[j] * r2 * si;
    pt[off++] = rli;
    pt[off++] = li;
    pt[off++] = 0;
    pt[off++] = 1 + dis;
    pt[off++] = 0.5;
    tmp = si * cd + co * sd;
    co = co * cd - si * sd;
    si = tmp;
    for (var j = 0; j < 3; j++)
      pt[off++] =
        b[j + 9] +
        lb * b[j] +
        bn[j + 3] +
        b[j] * (rm * co + 0.15) +
        bn[j] * rm * si;
    pt[off++] = rli;
    pt[off++] = li;
    pt[off++] = 0;
    pt[off++] = 0.75 + dis;
    pt[off++] = 0;
    tmp = si * cd + co * sd;
    co = co * cd - si * sd;
    si = tmp;
    for (var j = 0; j < 3; j++)
      pt[off++] =
        b[j + 9] +
        lb * b[j] +
        bn[j + 3] +
        b[j] * (r1 * co + 0.1) +
        bn[j] * r1 * si;
    pt[off++] = rli;
    pt[off++] = li;
    pt[off++] = 0;
    pt[off++] = 0.4 + dis;
    pt[off++] = 0;
  }
  for (var i = 0; i < 5; i++) {
    for (var j = 0; j < 4; j++) {
      ind[pi2++] = t0;
      ind[pi2++] = t++;
      ind[pi2++] = t;
    }
    t++;
  }
  dis = b[13];
  for (var j = 0; j < 3; j++) pt[off++] = b[j + 9] + R * b[j + 3];
  pt[off++] = 0.5;
  pt[off++] = 0.5;
  pt[off++] = 0;
  pt[off++] = dis;
  pt[off++] = 0;
  for (var j = 0; j < 3; j++) pt[off++] = b[j + 9] + R * b[j + 6];
  pt[off++] = 0.5;
  pt[off++] = 0.5;
  pt[off++] = 0;
  pt[off++] = dis;
  pt[off++] = 0;
  for (var j = 0; j < 3; j++) pt[off++] = b[j + 9] - R * b[j + 3];
  pt[off++] = 0.5;
  pt[off++] = 0.5;
  pt[off++] = 0;
  pt[off++] = dis;
  pt[off++] = 0;
  for (var j = 0; j < 3; j++) pt[off++] = b[j + 9] - R * b[j + 6];
  pt[off++] = 0.5;
  pt[off++] = 0.5;
  pt[off++] = 0;
  pt[off++] = dis;
  pt[off++] = 0;
  ind[pi2++] = t0;
  ind[pi2++] = t++;
  ind[pi2++] = t;
  ind[pi2++] = t0;
  ind[pi2++] = t++;
  ind[pi2++] = t;
  ind[pi2++] = t0;
  ind[pi2++] = t++;
  ind[pi2++] = t;
  ind[pi2++] = t0;
  ind[pi2++] = t;
  ind[pi2++] = t - 3;
}

function leaf3(b, R, lb, y0, ang, green, age) {
  console.log("age", age);
  var t0 = off / 8,
    r1 = 0.25 * age,
    rm = 0.39 * age,
    r2 = 0.5 * age,
    si = Math.sin(ang),
    co = Math.cos(ang),
    fi;
  var t = b[0] * co - b[2] * si;
  b[2] = b[0] * si + b[2] * co;
  b[0] = t;
  b[1] += y0;
  var bn = new Float32Array(6);
  var no = Math.sqrt(b[0] * b[0] + b[2] * b[2]);
  bn[0] = b[2] / no;
  bn[2] = -b[0] / no;
  bn[4] = 0.05 * 1.25; // change depth // previous val -0.05;
  for (var j = 0; j < 3; j++) pt[off++] = b[j + 9] + lb * b[j];
  var li = 0.02 * (1 - Math.abs(b[1])); // change shadow property of the leaf
  pt[off++] = 0;
  pt[off++] = 0.4 + li;
  pt[off++] = 0;
  pt[off++] = 0;
  pt[off++] = 0;

  t = t0 + 1;
  var sd = Math.sin(Math.PI * 0.3),
    cd = Math.cos(Math.PI * 0.3),
    si = 0,
    co = -1,
    tmp;
  for (var i = 0; i < 5; i++) {
    for (var j = 0; j < 3; j++)
      pt[off++] =
        b[j + 9] +
        lb * b[j] +
        bn[j + 3] +
        b[j] * (r1 * co + 0.1) +
        bn[j] * r1 * si;
    pt[off++] = 0;
    pt[off++] = 0.9 + li;
    pt[off++] = 0;
    pt[off++] = 0.125;
    pt[off++] = 1;
    tmp = si * cd + co * sd;
    co = co * cd - si * sd;
    si = tmp;
    for (var j = 0; j < 3; j++)
      pt[off++] =
        b[j + 9] +
        lb * b[j] +
        bn[j + 3] +
        b[j] * (rm * co + 0.15) +
        bn[j] * rm * si;
    pt[off++] = 0;
    pt[off++] = 0.6 + li;
    pt[off++] = 0;
    pt[off++] = 0.5;
    pt[off++] = 1;
    tmp = si * cd + co * sd;
    co = co * cd - si * sd;
    si = tmp;
    for (var j = 0; j < 3; j++)
      pt[off++] =
        b[j + 9] +
        lb * b[j] +
        bn[j + 3] +
        b[j] * (r2 * co + 0.2) +
        bn[j] * r2 * si;
    pt[off++] = 0;
    pt[off++] = 0.4 + li;
    pt[off++] = 0;
    pt[off++] = 1;
    pt[off++] = 1;
    tmp = si * cd + co * sd;
    co = co * cd - si * sd;
    si = tmp;
    for (var j = 0; j < 3; j++)
      pt[off++] =
        b[j + 9] +
        lb * b[j] +
        bn[j + 3] +
        b[j] * (rm * co + 0.15) +
        bn[j] * rm * si;
    pt[off++] = 0;
    pt[off++] = 0.6 + li;
    pt[off++] = 0;
    pt[off++] = 1;
    pt[off++] = 0.5;
    tmp = si * cd + co * sd;
    co = co * cd - si * sd;
    si = tmp;
    for (var j = 0; j < 3; j++)
      pt[off++] =
        b[j + 9] +
        lb * b[j] +
        bn[j + 3] +
        b[j] * (r1 * co + 0.1) +
        bn[j] * r1 * si;
    pt[off++] = 0;
    pt[off++] = 0.9 + li; // innard colors
    pt[off++] = 0;
    pt[off++] = 1;
    pt[off++] = 0.125;
  }
  for (var i = 0; i < 5; i++) {
    for (var j = 0; j < 4; j++) {
      ind[pi2++] = t0;
      ind[pi2++] = t++;
      ind[pi2++] = t;
    }
    t++;
  }
  for (var j = 0; j < 3; j++) pt[off++] = b[j + 9] + R * b[j + 3];
  pt[off++] = 0.5;
  pt[off++] = 0.5;
  pt[off++] = pt[off++] = pt[off++] = 0;
  for (var j = 0; j < 3; j++) pt[off++] = b[j + 9] + R * b[j + 6];
  pt[off++] = 0.5;
  pt[off++] = 0.5;
  pt[off++] = pt[off++] = pt[off++] = 0;
  for (var j = 0; j < 3; j++) pt[off++] = b[j + 9] - R * b[j + 3];
  pt[off++] = 0.5;
  pt[off++] = 0.5;
  pt[off++] = pt[off++] = pt[off++] = 0;
  for (var j = 0; j < 3; j++) pt[off++] = b[j + 9] - R * b[j + 6];
  pt[off++] = 0.5;
  pt[off++] = 0.5;
  pt[off++] = pt[off++] = pt[off++] = 0;
  ind[pi2++] = t0;
  ind[pi2++] = t++;
  ind[pi2++] = t;
  ind[pi2++] = t0;
  ind[pi2++] = t++;
  ind[pi2++] = t;
  ind[pi2++] = t0;
  ind[pi2++] = t++;
  ind[pi2++] = t;
  ind[pi2++] = t0;
  ind[pi2++] = t;
  ind[pi2++] = t - 3;
}

function twig(b, th, nfi, R0, R2, l0, l1, stu, sub) {
  var t = off / 8,
    ssu = (su - 1) / stu + 1,
    tm = t + 2 * ssu - 1;
  for (var j = 0; j < sv - 1; j++) {
    for (var u = 0; u < 2 * ssu - 2; u++) {
      ind[pi++] = t++;
      ind[pi++] = t;
      ind[pi++] = tm;
      ind[pi++] = tm++;
      ind[pi++] = t;
      ind[pi++] = tm;
    }
    t++;
    tm++;
  }
  var co = Math.cos(th),
    si = Math.sin(th),
    ta2 = Math.tan(th * 0.5),
    fi = (Math.PI * nfi) / (ssu - 1),
    cof = Math.cos(fi),
    sif = Math.sin(fi);
  var cp = new Float32Array(6),
    bn = new Float32Array(6);
  for (var i = 0; i < 3; i++) {
    cp[i] = b[i + 9];
    cp[i + 3] = b[i + 9] + l0 * b[i];
    bn[i] = sif * b[i + 3] + cof * b[i + 6]; // fi rot
    b[i + 3] = cof * b[i + 3] - sif * b[i + 6];
    bn[i + 3] = ta2 * b[i] + bn[i];
    b[i + 6] = si * b[i] + co * bn[i]; // th rot
    b[i] = co * b[i] - si * bn[i];
    b[i + 9] = cp[i + 3] + l1 * b[i];
  }
  b[12] -= nfi / (ssu - 1);
  var ls = 0.5 * (l0 + l1),
    R1 = (l0 * R0 + l1 * R2) / (l0 + l1);
  for (var v = 0; v < sv; v++) {
    var t = v / (sv - 1),
      t1 = 1 - t,
      B0 = t1 * t1,
      B1 = 2 * t * t1,
      B2 = t * t;
    for (var u = 0, u1 = 0; u < 2 * su - 1; u += stu, u1++) {
      var si = sin[u],
        co = cos[u],
        si2 = si,
        co2 = co;
      if (sub && u1 & 1) {
        si2 = (sin[u - stu] + sin[u + stu]) * 0.5;
        co2 = (cos[u - stu] + cos[u + stu]) * 0.5;
      }
      for (var i = 0; i < 3; i++) {
        pt[off++] =
          (cp[i] + (bn[i] * si + b[i + 3] * co) * R0) * B0 +
          (cp[i + 3] + (bn[i + 3] * si + b[i + 3] * co) * R1) * B1 +
          (b[i + 9] + (b[i + 6] * si2 + b[i + 3] * co2) * R2) * B2;
      }
      var li = 0.5 * (b[3] * co + (B0 * bn[0] + B1 * bn[3] + B2 * b[6]) * si);
      if (li < 0) li = 0;
      //      pt[off++] = .7 - .3*(u1 & 1); pt[off++] = .6 - .3*(u1 & 1); pt[off++] = 0;
      pt[off++] = 0.6 + li;
      pt[off++] = 0.3 + li;
      pt[off++] = 0;
      pt[off++] = ls * t + b[13];
      pt[off++] = u / (su - 1) + b[12];
    }
  }
  b[13] += ls;
}
function bark_tex(k, img) {
  var kk = k * k,
    k1 = k - 1,
    kk1 = kk - 1;
  var hk = new Float32Array(k * k);
  var st = k / 4,
    z0 = 1;
  while (st > 1) {
    var st2 = st >> 1,
      stk = st << 8,
      stk2 = stk >> 1;
    for (var j = 0; j < kk; j += stk) {
      var jp = (j + stk) & kk1;
      for (var i = 0; i < k; i += st) {
        var ip = (i + st) & k1;
        var h0 = hk[i + j],
          h1 = hk[ip + j],
          h2 = hk[i + jp];
        hk[i + st2 + j + stk2] =
          (h0 + h1 + h2 + hk[ip + jp]) * 0.25 + z0 * (Math.random() - 0.5);
        hk[i + st2 + j] = (h0 + h1) * 0.5 + z0 * (Math.random() - 0.5);
        hk[i + j + stk2] = (h0 + h2) * 0.5 + z0 * (Math.random() - 0.5);
      }
    }
    st = st2;
    z0 *= 0.55;
  }
  var t = 0,
    s = 0;
  for (var i = 0; i < k; i++)
    for (var j = 0; j < k; j++) {
      img[t++] = 800 * Math.abs(hk[s++] + 0.6 * Math.sin((Math.PI * i) / 12)); // change bark tex
      img[t++] = img[t++] = img[t++] = 0;
    }
}
function leaf_tex(k, img) {
  var kk = k * k,
    k1 = k - 1,
    a = 0.7 / k;
  var hk = new Float32Array(k * k);
  for (var j = 1; j < k - 1; j++) hk[j * (k + 1)] = 1 - a * j;
  for (var j = 1; j < k / 2 - 1; j++) {
    hk[kk / 2 - j] = a * j;
    hk[kk / 2 + k / 2 + j * k] = 0.5 - a * j;
  }
  for (var j = 1; j < (3 * k) / 4 - 1; j++) {
    hk[kk / 4 - j] = a * j;
    hk[kk / 4 + k / 4 + j * k] = 0.75 - a * j;
  }
  for (var it = 0; it < 5; it++)
    for (var j = k; j < k * k1; j++)
      hk[j] =
        (4 * hk[j] + hk[j + 1] + hk[j - 1] + hk[j + k] + hk[j - k]) * 0.125;
  var t = 0;
  for (var i = 0; i < kk; i++) {
    img[t++] = 255 * (1 - hk[i]);
    img[t++] = img[t++] = img[t++] = 0;
  }
}
function leaf_tex3(k, img) {
  var kk = k * k,
    k2 = k / 2,
    k4 = k / 4,
    kk2 = k2 * k,
    a = 1 / k;
  var hk = new Float32Array(kk);
  for (var j = 1; j < k - 1; j++) hk[kk2 + j] = 1 - a * j;
  a *= 4;
  for (var j = 1; j < k / 4 - 1; j++)
    hk[kk2 + k2 - j * (k + k - 1)] =
      hk[kk2 + k2 + j * (k + k + 1)] =
      hk[kk2 + k4 + j * (k + k + 1)] =
      hk[kk2 + k4 - j * (k + k - 1)] =
        1 - a * j;
  for (var it = 0; it < 5; it++)
    for (var j = k; j < k * (k - 1); j++)
      hk[j] =
        (4 * hk[j] + hk[j + 1] + hk[j - 1] + hk[j + k] + hk[j - k]) * 0.125;
  var t = 0;
  for (var i = 0; i < kk; i++) {
    img[t++] = 255 * (1 - hk[i]);
    img[t++] = img[t++] = img[t++] = 0;
  }
}

/*
  function leaves2(b, R) {
    if (b[1] > 0) {
      b[1] = b[1] * (1 - b[1] * b[1]);
      var no = 0;
      for (var j = 0; j < 3; j++) no += b[j] * b[j];
      no = 1 / Math.sqrt(no);
      for (var j = 0; j < 3; j++) b[j] *= no;
    }
    leaf2(b, R, 0.8, 0.1, 0);
    leaf2(b, R, 0.7, -0.1, 0.9);
    leaf2(b, R, 0.7, 0.1, -1.8);
    leaf2(b, R, 1.4, -0.15, 0.4);
    leaf2(b, R, 1.5, -0.05, 0.9);
  }
  function leaf2(b, R, lb, y0, ang) {
    var t0 = off / 8,
      r1 = 0.25,
      rm = 0.39,
      r2 = 0.5,
      si = Math.sin(ang),
      co = Math.cos(ang),
      fi;
    var t = b[0] * co - b[2] * si;
    b[2] = b[0] * si + b[2] * co;
    b[0] = t;
    b[1] += y0;
    var bn = new Float32Array(6);
    var no = Math.sqrt(b[0] * b[0] + b[2] * b[2]);
    bn[0] = b[2] / no;
    bn[2] = -b[0] / no;
    bn[4] = -0.05;
    for (var j = 0; j < 3; j++) pt[off++] = b[j + 9] + lb * b[j];
    var li = 0.2 * (1 - Math.abs(b[1]));
    pt[off++] = 0;
    pt[off++] = 0.3 + li;
    pt[off++] = 0;
    pt[off++] = 0;
    pt[off++] = 0;
    t = t0 + 1;
    var sd = Math.sin(Math.PI * 0.08),
      cd = Math.cos(Math.PI * 0.08),
      si = si / 2,
      co = -1,
      tmp;
    for (var i = 0; i < 5; i++) {
      tmp = si * cd + co * sd;
      co = co * cd - si * sd;
      si = tmp;
      for (var j = 0; j < 3; j++)
        pt[off++] =
          b[j + 9] +
          lb * b[j] +
          bn[j + 3] +
          b[j] * (r1 * co + 0.1) +
          bn[j] * r1 * si;
      pt[off++] = 0;
      pt[off++] = 0.7 + li;
      pt[off++] = 0;
      pt[off++] = 0.25;
      pt[off++] = 1;
      tmp = si * cd + co * sd;
      co = co * cd - si * sd;
      si = tmp;
      for (var j = 0; j < 3; j++)
        pt[off++] =
          b[j + 9] +
          lb * b[j] +
          bn[j + 3] +
          b[j] * (rm * co + 0.15) +
          bn[j] * rm * si;
      pt[off++] = 0;
      pt[off++] = 0.7 + li;
      pt[off++] = 0;
      pt[off++] = 0.5;
      pt[off++] = 1;
      tmp = si * cd + co * sd;
      co = co * cd - si * sd;
      si = tmp;
      for (var j = 0; j < 3; j++)
        pt[off++] =
          b[j + 9] +
          lb * b[j] +
          bn[j + 3] +
          b[j] * (r2 * co + 0.2) +
          bn[j] * r2 * si;
      pt[off++] = 0;
      pt[off++] = 0.6 + li;
      pt[off++] = 0;
      pt[off++] = 1;
      pt[off++] = 1;
      tmp = si * cd + co * sd;
      co = co * cd - si * sd;
      si = tmp;
      for (var j = 0; j < 3; j++)
        pt[off++] =
          b[j + 9] +
          lb * b[j] +
          bn[j + 3] +
          b[j] * (rm * co + 0.15) +
          bn[j] * rm * si;
      pt[off++] = 0;
      pt[off++] = 0.7 + li;
      pt[off++] = 0;
      pt[off++] = 1;
      pt[off++] = 0.5;
      tmp = si * cd + co * sd;
      co = co * cd - si * sd;
      si = tmp;
      for (var j = 0; j < 3; j++)
        pt[off++] =
          b[j + 9] +
          lb * b[j] +
          bn[j + 3] +
          b[j] * (r1 * co + 0.1) +
          bn[j] * r1 * si;
      pt[off++] = 0;
      pt[off++] = 0.7 + li;
      pt[off++] = 0;
      pt[off++] = 1;
      pt[off++] = 0.25;
    }
    for (var i = 0; i < 5; i++) {
      for (var j = 0; j < 4; j++) {
        ind[pi2++] = t0;
        ind[pi2++] = t++;
        ind[pi2++] = t;
      }
      t++;
    }
    for (var j = 0; j < 3; j++) pt[off++] = b[j + 9] + R * b[j + 3];
    pt[off++] = 0.5;
    pt[off++] = 0.5;
    pt[off++] = pt[off++] = pt[off++] = 0;
    for (var j = 0; j < 3; j++) pt[off++] = b[j + 9] + R * b[j + 6];
    pt[off++] = 0.5;
    pt[off++] = 0.5;
    pt[off++] = pt[off++] = pt[off++] = 0;
    for (var j = 0; j < 3; j++) pt[off++] = b[j + 9] - R * b[j + 3];
    pt[off++] = 0.5;
    pt[off++] = 0.5;
    pt[off++] = pt[off++] = pt[off++] = 0;
    for (var j = 0; j < 3; j++) pt[off++] = b[j + 9] - R * b[j + 6];
    pt[off++] = 0.5;
    pt[off++] = 0.5;
    pt[off++] = pt[off++] = pt[off++] = 0;
    ind[pi2++] = t0;
    ind[pi2++] = t++;
    ind[pi2++] = t;
    ind[pi2++] = t0;
    ind[pi2++] = t++;
    ind[pi2++] = t;
    ind[pi2++] = t0;
    ind[pi2++] = t++;
    ind[pi2++] = t;
    ind[pi2++] = t0;
    ind[pi2++] = t;
    ind[pi2++] = t - 3;
  }
  */
