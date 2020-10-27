//==========================================================================//
// TFER_PMA CODE -----------------------------------------------------------//
var pi = 3.14159265359

var prop_pma = function(opts) {
  if (typeof opts == 'undefined' || (opts != null && opts.hasOwnProperty("__kwargtrans__"))) {
    ;
    var opts = 'Olfert';
  };
  if (opts == 'Olfert') {
    var prop = {
      'r1': 0.06,
      'r2': 0.061,
      'L': 0.2,
      'p': 1,
      'T': 293,
      'Q': (3 / 1000) / 60,
      'omega_hat': 32 / 33
    };
  } else if (opts == 'Buckley') {
    var prop = {
      'r1': 0.025,
      'r2': 0.024,
      'L': 0.1,
      'omega': ((13350 * 2) * pi) / 60,
      'p': 1,
      'T': 298,
      'Q': 0.00102 / 60,
      'omega_hat': 1
    };
  }
  prop['rc'] = (prop['r1'] + prop['r2']) / 2;
  prop['r_hat'] = prop['r1'] / prop['r2'];
  prop['del'] = (prop['r2'] - prop['r1']) / 2;
  prop['A'] = pi * (Math.pow(prop['r2'], 2) - Math.pow(prop['r1'], 2));
  prop['v_bar'] = prop['Q'] / prop['A'];
  var kB = 1.3806488e-23;
  prop['D'] = (function __lambda__(B) {
    return (kB * prop['T']) * B;
  });
  prop['Dm'] = 3;
  prop['m0'] = 4.7124e-25; // mass @ dm = 1 in fg
  return prop;
};


var get_setpoint = function(prop) {
  if (typeof prop == 'undefined' || (prop != null && prop.hasOwnProperty("__kwargtrans__"))) {
    ;
    var prop = {};
  };
  if (!(prop)) {
    var prop = prop_pma();
  }
  var sp = {
    'm_star': null,
    'V': null,
    'omega': null,
    'omega1': null,
    'Rm': null
  };
  if (arguments.length == 3) {
    sp[arguments[1]] = arguments[2];
    sp['Rm'] = 3;
  } else if (arguments.length == 5) {
    sp[arguments[1]] = arguments[2];
    sp[arguments[3]] = arguments[4];
  }
  var e = 1.60218e-19;
  if (sp['m_star'] == null) {
    if (!(sp['omega1'])) {
      sp['omega1'] = sp['omega'] / ((Math.pow(prop['r_hat'], 2) - prop['omega_hat']) / (Math.pow(prop['r_hat'], 2) - 1) + ((Math.pow(prop['r1'], 2) * (prop['omega_hat'] - 1)) / (Math.pow(prop['r_hat'], 2) - 1)) / Math.pow(prop['rc'], 2));
    }
    sp['alpha'] = (sp['omega1'] * (Math.pow(prop['r_hat'], 2) - prop['omega_hat'])) / (Math.pow(prop['r_hat'], 2) - 1);
    sp['beta'] = ((sp['omega1'] * Math.pow(prop['r1'], 2)) * (prop['omega_hat'] - 1)) / (Math.pow(prop['r_hat'], 2) - 1);
    sp['m_star'] = sp['V'] / ((Math.log(1 / prop['r_hat']) / e) * Math.pow(sp['alpha'] * prop['rc'] + sp['beta'] / prop['rc'], 2));
    sp['omega'] = sp['alpha'] + sp['beta'] / Math.pow(prop['rc'], 2);
    sp['omega2'] = sp['alpha'] + sp['beta'] / Math.pow(prop['r2'], 2);
  } else if (sp['omega1'] != null) {
    sp['alpha'] = (sp['omega1'] * (Math.pow(prop['r_hat'], 2) - prop['omega_hat'])) / (Math.pow(prop['r_hat'], 2) - 1);
    sp['beta'] = ((sp['omega1'] * Math.pow(prop['r1'], 2)) * (prop['omega_hat'] - 1)) / (Math.pow(prop['r_hat'], 2) - 1);
    sp['V'] = ((sp['m_star'] * Math.log(1 / prop['r_hat'])) / e) * Math.pow(sp['alpha'] * prop['rc'] + sp['beta'] / prop['rc'], 2);
    sp['omega2'] = sp['alpha'] + sp['beta'] / Math.pow(prop['r2'], 2);
    sp['omega'] = sp['alpha'] + sp['beta'] / Math.pow(prop['rc'], 2);
  } else if (sp['omega'] != null) {
    sp['omega1'] = sp['omega'] / ((Math.pow(prop['r_hat'], 2) - prop['omega_hat']) / (Math.pow(prop['r_hat'], 2) - 1) + ((Math.pow(prop['r1'], 2) * (prop['omega_hat'] - 1)) / (Math.pow(prop['r_hat'], 2) - 1)) / Math.pow(prop['rc'], 2));
    sp['alpha'] = (sp['omega1'] * (Math.pow(prop['r_hat'], 2) - prop['omega_hat'])) / (Math.pow(prop['r_hat'], 2) - 1);
    sp['beta'] = ((sp['omega1'] * Math.pow(prop['r1'], 2)) * (prop['omega_hat'] - 1)) / (Math.pow(prop['r_hat'], 2) - 1);
    sp['V'] = ((sp['m_star'] * Math.log(1 / prop['r_hat'])) / e) * Math.pow(sp['alpha'] * prop['rc'] + sp['beta'] / prop['rc'], 2);
    sp['omega2'] = sp['alpha'] + sp['beta'] / Math.pow(prop['r2'], 2);
  } else if (sp['V'] != null) {
    var v_theta_rc = Math.sqrt((sp['V'] * e) / (sp['m_star'] * Math.log(1 / prop['r_hat'])));
    var A = (prop['rc'] * (Math.pow(prop['r_hat'], 2) - prop['omega_hat'])) / (Math.pow(prop['r_hat'], 2) - 1) + (1 / prop['rc']) * ((Math.pow(prop['r1'], 2) * (prop['omega_hat'] - 1)) / (Math.pow(prop['r_hat'], 2) - 1));
    sp['omega1'] = v_theta_rc / A;
    sp['alpha'] = (sp['omega1'] * (Math.pow(prop['r_hat'], 2) - prop['omega_hat'])) / (Math.pow(prop['r_hat'], 2) - 1);
    sp['beta'] = ((sp['omega1'] * Math.pow(prop['r1'], 2)) * (prop['omega_hat'] - 1)) / (Math.pow(prop['r_hat'], 2) - 1);
    sp['omega2'] = sp['alpha'] + sp['beta'] / Math.pow(prop['r2'], 2);
    sp['omega'] = sp['alpha'] + sp['beta'] / Math.pow(prop['rc'], 2);
  } else if (sp['Rm'] != null) {
    var n_B = get_nb(sp['m_star'], prop);
    var __left0__ = mp2zp(sp['m_star'], 1, prop['T'], prop['p'], prop);
    var B_star = __left0__[0];
    sp['m_max'] = sp['m_star'] * (1 / sp['Rm'] + 1);
    sp['omega'] = Math.sqrt(prop['Q'] / ((((((sp['m_star'] * B_star) * 2) * pi) * Math.pow(prop['rc'], 2)) * prop['L']) * (Math.pow(sp['m_max'] / sp['m_star'], n_B + 1) - Math.pow(sp['m_max'] / sp['m_star'], n_B))));
    sp['omega1'] = sp['omega'] / ((Math.pow(prop['r_hat'], 2) - prop['omega_hat']) / (Math.pow(prop['r_hat'], 2) - 1) + ((Math.pow(prop['r1'], 2) * (prop['omega_hat'] - 1)) / (Math.pow(prop['r_hat'], 2) - 1)) / Math.pow(prop['rc'], 2));
    sp['alpha'] = (sp['omega1'] * (Math.pow(prop['r_hat'], 2) - prop['omega_hat'])) / (Math.pow(prop['r_hat'], 2) - 1);
    sp['beta'] = ((sp['omega1'] * Math.pow(prop['r1'], 2)) * (prop['omega_hat'] - 1)) / (Math.pow(prop['r_hat'], 2) - 1);
    sp['omega2'] = sp['alpha'] + sp['beta'] / Math.pow(prop['r2'], 2);
    sp['V'] = ((sp['m_star'] * Math.log(1 / prop['r_hat'])) / e) * Math.pow(sp['alpha'] * prop['rc'] + sp['beta'] / prop['rc'], 2);
  } else {
    console.log('No setpoint specified.');
  }
  var m_star = sp['m_star'];
  if (!(sp['Rm'])) {
    var __left0__ = get_resolution(sp['m_star'], sp['omega'], prop);
    sp['Rm'] = __left0__[0];
    sp['m_max'] = __left0__[1];
  }
  return [sp, m_star];
};

var get_resolution = function(m_star, omega, prop) {
  print('Finding resolution...');
  var n_B = get_nb(m_star, prop);
  var __left0__ = mp2zp(m_star, 1, prop['T'], prop['p'], prop);
  var B_star = __left0__[0];
  var t0 = prop['Q'] / ((((((m_star * B_star) * 2) * pi) * prop['L']) * Math.pow(omega, 2)) * Math.pow(prop['rc'], 2));
  var m_rat = (function __lambda__(Rm) {
    return 1 / Rm + 1;
  });
  var fun = (function __lambda__(Rm) {
    return Math.pow(m_rat(Rm), n_B + 1) - Math.pow(m_rat(Rm), n_B);
  });
  var Rm = scipy.optimize.fmin((function __lambda__(Rm) {
    return Math.pow(t0 - fun(Rm), 2);
  }), __kwargtrans__({
    x0: 5
  }));
  var Rm = Rm[0];
  var m_max = m_star * (1 / Rm + 1);
  print('Complete.');
  print(' ');
  return [Rm, m_max];
};

var get_nb = function(m_star, prop) {
  var m_high = m_star * 1.001;
  var m_low = m_star * 0.999;
  var __left0__ = mp2zp(m_high, 1, prop['T'], prop['p'], prop);
  var B_high = __left0__[0];
  var __left0__ = mp2zp(m_low, 1, prop['T'], prop['p'], prop);
  var B_low = __left0__[0];
  var n_B = Math.log10(B_high / B_low) / Math.log10(m_high / m_low);
  return n_B;
};

var dm2zp = function(d, z, T, p) {
  if (typeof T == 'undefined' || (T != null && T.hasOwnProperty("__kwargtrans__"))) {
    ;
    var T = 0.0;
  };
  if (typeof p == 'undefined' || (p != null && p.hasOwnProperty("__kwargtrans__"))) {
    ;
    var p = 0.0;
  };
  var e = 1.6022e-19;
  if (T == 0.0 || p == 0.0) {
    var mu = 1.82e-05;
    var B = Cc(d) / (((3 * pi) * mu) * d);
  } else {
    var S = 110.4;
    var T_0 = 296.15;
    var vis_23 = 1.83245e-05;
    var mu = (vis_23 * Math.pow(T / T_0, 1.5)) * ((T_0 + S) / (T + S));
    var B = Cc(d, T, p) / (((3 * pi) * mu) * d);
  }
  var Zp = (B * e) * z;
  return [B, Zp];
};

var Cc = function(d, T, p) {
  if (typeof T == 'undefined' || (T != null && T.hasOwnProperty("__kwargtrans__"))) {
    ;
    var T = 0.0;
  };
  if (typeof p == 'undefined' || (p != null && p.hasOwnProperty("__kwargtrans__"))) {
    ;
    var p = 0.0;
  };
  if (T == 0.0 || p == 0.0) {
    var mfp = 6.65e-08;
    var A1 = 1.257;
    var A2 = 0.4;
    var A3 = 0.55;
  } else {
    var S = 110.4;
    var mfp_0 = 6.73e-08;
    var T_0 = 296.15;
    var p_0 = 101325;
    var p = p * p_0;
    var mfp = ((mfp_0 * Math.pow(T / T_0, 2)) * (p_0 / p)) * ((T_0 + S) / (T + S));
    var A1 = 1.165;
    var A2 = 0.483;
    var A3 = 0.997 / 2;
  }
  var Kn = (2 * mfp) / d;
  var Cc = 1 + Kn * (A1 + A2 * Math.exp(-(2 * A3) / Kn));
  return Cc;
};

var mp2zp = function(m, z, T, P, prop) {
  if (typeof T == 'undefined' || (T != null && T.hasOwnProperty("__kwargtrans__"))) {
    ;
    var T = 0.0;
  };
  if (typeof P == 'undefined' || (P != null && P.hasOwnProperty("__kwargtrans__"))) {
    ;
    var P = 0.0;
  };
  if (typeof prop == 'undefined' || (prop != null && prop.hasOwnProperty("__kwargtrans__"))) {
    ;
    var prop = {};
  };
  if (!('m0' in prop) || !('Dm' in prop)) {
    sys.exit('Please specify the mass-mobility relation parameters in prop.');
  }
  var d = Math.pow(m / prop['m0'], 1 / prop['Dm']) * 1e-9;
  if (T == 0.0 || P == 0.0) {
    var __left0__ = dm2zp(d, z);
    var B = __left0__[0];
    var Zp = __left0__[1];
  } else {
    var __left0__ = dm2zp(d, z, T, P);
    var B = __left0__[0];
    var Zp = __left0__[1];
  }
  return [B, Zp, d];
};


var tfer_1S = function(sp, m, d, z, prop) {
  if (typeof prop == 'undefined' || (prop != null && prop.hasOwnProperty("__kwargtrans__"))) {
    ;
    var prop = {};
  };
  var __left0__ = parse_inputs(sp, m, d, z, prop);
  var tau = __left0__[0];
  var rs = __left0__[3];
  var lam = (((2 * tau) * (Math.pow(sp['alpha'], 2) - Math.pow(sp['beta'], 2) / Math.pow(rs, 4))) * prop['L']) / prop['v_bar'];
  var G0 = (function __lambda__(r) {
    return rs + (r - rs) * Math.exp(-(lam));
  });
  var ra = Math.min(prop['r2'], Math.max(prop['r1'], G0(prop['r1'])));
  var rb = Math.min(prop['r2'], Math.max(prop['r1'], G0(prop['r2'])));
  var Lambda = (1 / (2 * prop['del'])) * (rb - ra);
  return [Lambda, G0];
};

var tfer_1C = function(sp, m, d, z, prop) {
  if (typeof prop == 'undefined' || (prop != null && prop.hasOwnProperty("__kwargtrans__"))) {
    ;
    var prop = {};
  };
  var __left0__ = parse_inputs(sp, m, d, z, prop);
  var tau = __left0__[0];
  var C0 = __left0__[1];
  var C3 = tau * (((Math.pow(sp['alpha'], 2) * prop['rc'] + ((2 * sp['alpha']) * sp['beta']) / prop['rc']) + Math.pow(sp['beta'], 2) / Math.pow(prop['rc'], 3)) - C0 / (m * prop['rc']));
  var C4 = tau * (((Math.pow(sp['alpha'], 2) - ((2 * sp['alpha']) * sp['beta']) / Math.pow(prop['rc'], 2)) - (3 * Math.pow(sp['beta'], 2)) / Math.pow(prop['rc'], 4)) + C0 / (m * Math.pow(prop['rc'], 2)));
  var G0 = (function __lambda__(r) {
    return (prop['rc'] + ((r - prop['rc']) + C3 / C4) * Math.exp((-(C4) * prop['L']) / prop['v_bar'])) - C3 / C4;
  });
  var ra = Math.min(prop['r2'], Math.max(prop['r1'], G0(prop['r1'])));
  var rb = Math.min(prop['r2'], Math.max(prop['r1'], G0(prop['r2'])));
  var Lambda = (1 / (2 * prop['del'])) * (rb - ra);
  return [Lambda, G0];
};

var tfer_1C_diff = function(sp, m, d, z, prop) {
  if (typeof prop == 'undefined' || (prop != null && prop.hasOwnProperty("__kwargtrans__"))) {
    ;
    var prop = {};
  };
  var __left0__ = parse_inputs(sp, m, d, z, prop);
  var D = __left0__[2];
  var sig = Math.sqrt(((2 * prop['L']) * D) / prop['v_bar']);
  var __left0__ = tfer_1C(sp, m, d, z, prop);
  var G0 = __left0__[1];
  var rho_fun = (function __lambda__(G, r) {
    return (G - r) / (Math.sqrt(2) * sig);
  });
  var kap_fun = (function __lambda__(G, r) {
    return (G - r) * math.erf(rho_fun(G, r)) + (sig * Math.sqrt(2 / pi)) * Math.exp(-(Math.pow(rho_fun(G, r), 2)));
  });
  var K22 = kap_fun(G0(prop['r2']), prop['r2']);
  var K21 = kap_fun(G0(prop['r2']), prop['r1']);
  var K12 = kap_fun(G0(prop['r1']), prop['r2']);
  var K11 = kap_fun(G0(prop['r1']), prop['r1']);
  var Lambda = (-(1) / (4 * prop['del'])) * (((K22 - K12) - K21) + K11);
  Lambda[K22 > 100.0] = 0;
  Lambda[Math.abs(Lambda) < 1e-10] = 0;
  return [Lambda, G0];
};

var tfer_W1 = function(sp, m, d, z, prop) {
  var __left0__ = parse_inputs(sp, m, d, z, prop);
  var tau = __left0__[0];
  var C0 = __left0__[1];
  var rs = __left0__[3];
  var lam = (((2 * tau) * (Math.pow(sp['alpha'], 2) - Math.pow(sp['beta'], 2) / Math.pow(rs, 4))) * prop['L']) / prop['v_bar'];
  var G0 = (function __lambda__(r) {
    return (1 / (sp['omega1'] * Math.sqrt(m))) * Math.sqrt(((m * Math.pow(sp['omega1'], 2)) * Math.pow(r, 2) - C0) * Math.exp(-(lam)) + C0);
  });
  var ra = Math.min(prop['r2'], Math.max(prop['r1'], G0(prop['r1'])));
  var rb = Math.min(prop['r2'], Math.max(prop['r1'], G0(prop['r2'])));
  var Lambda = (1 / (2 * prop['del'])) * (rb - ra);
  return [Lambda, G0];
};

var tfer_W1_diff = function(sp, m, d, z, prop) {
  var __left0__ = parse_inputs(sp, m, d, z, prop);
  var D = __left0__[2];
  var sig = Math.sqrt(((2 * prop['L']) * D) / prop['v_bar']);
  var __left0__ = tfer_W1(sp, m, d, z, prop);
  var G0 = __left0__[1];
  var rho_fun = (function __lambda__(G, r) {
    return (G - r) / (Math.sqrt(2) * sig);
  });
  var kap_fun = (function __lambda__(G, r) {
    return (G - r) * math.erf(rho_fun(G, r)) + (sig * Math.sqrt(2 / pi)) * Math.exp(-(Math.pow(rho_fun(G, r), 2)));
  });
  var K22 = kap_fun(G0(prop['r2']), prop['r2']);
  var K21 = kap_fun(G0(prop['r2']), prop['r1']);
  var K12 = kap_fun(G0(prop['r1']), prop['r2']);
  var K11 = kap_fun(G0(prop['r1']), prop['r1']);
  var Lambda = (-(1) / (4 * prop['del'])) * (((K22 - K12) - K21) + K11);
  Lambda[K22 > 100.0] = 0;
  Lambda[Math.abs(Lambda) < 1e-10] = 0;
  return [Lambda, G0];
};

var parse_inputs = function(sp, m, d, z, prop) {
  if (typeof d == 'undefined' || (d != null && d.hasOwnProperty("__kwargtrans__"))) {
    ;
    var d = 0;
  };
  if (typeof z == 'undefined' || (z != null && z.hasOwnProperty("__kwargtrans__"))) {
    ;
    var z = 1;
  };
  if (typeof prop == 'undefined' || (prop != null && prop.hasOwnProperty("__kwargtrans__"))) {
    ;
    var prop = {};
  };
  var e = 1.60218e-19;
  var q = z * e;
  if (d == 0 || d == null) {
    console.log('Invoking mass-mobility relation to determine Zp.');
    var __left0__ = mp2zp(m, z, prop['T'], prop['p'], prop);
    var B = __left0__[0];
  } else {
    var __left0__ = dm2zp(d, z, prop['T'], prop['p']);
    var B = __left0__[0];
  }
  var tau = B * m;
  var D = prop['D'](B) * z;
  var C0 = (sp['V'] * q) / Math.log(1 / prop['r_hat']);

  // Note: Whether to pick the +ive of -ive root for rs is chosen based on a
  // heuristic approach. Specifically, the root closer to the centerline is
  // chosen, except when the -ive root is zero (which is the case for APM
  // conditions, where the +ive root should always be used).
  // evaluate +ive and -ive roots
  r_m = (math.sqrt(C0 / m) - math.sqrt(C0 / m - 4 * sp['alpha'] * sp['beta'])) / (2 * sp['alpha']);
  r_p = (math.sqrt(C0 / m) + math.sqrt(C0 / m - 4 * sp['alpha'] * sp['beta'])) / (2 * sp['alpha']);

  // determine which root is closer to centerline radius
  function indexOfSmallest(a) {
    var lowest = 0;
    for (var i = 1; i < a.length; i++) {
      if (a[i] < a[lowest]) lowest = i;
    }
    return lowest;
  }
  idx = indexOfSmallest([Math.abs(r_m - prop['rc']), Math.abs(r_p - prop['rc'])]);
  if (r_m == 0) {
    idx = 2; // avoid zero values for APM case
  }

  // assign one of the two roots to rs
  rs = r_m; // by default use -ive root
  if (idx == 2) {
    rs = r_p // if closer to +ive root, use +ive root
  }

  // zero out cases where no equilibrium radius exists (also removes complex numbers)
  if (C0 / m < (4 * sp['alpha'] * sp['beta'])) {
    rs = 0;
  }

  rs = Math.max(rs, 1e-20) // avoids division by zero for z = 0 case

  return [tau, C0, D, rs];
};

var linspace = function(a, b, n) {
  if (typeof n === "undefined") n = Math.max(Math.round(b - a) + 1, 1);
  if (n < 2) {
    return n === 1 ? [a] : [];
  }
  var i, ret = Array(n);
  n--;
  for (i = n; i >= 0; i--) {
    ret[i] = (i * b + (n - i) * a) / n;
  }
  return ret;
};

var parse_fun = function(sp, m, d, prop, fun) {
  var Lambda = Array(m.length)
  for (ii in m) { // loop over particle mass
    Lambda[ii] = 0 // initialize at zero
    for (zz in z_vec) { // loop over integer charge states
      __left0__ = fun(sp, m[ii], d[ii], z_vec[zz], prop)
      Lambda[ii] = Lambda[ii] + __left0__[0]
    }
  }
  return Lambda;
};

console.log('Load complete.')
console.log(' ')

// END TFER_PMA CODE -------------------------------------------------------//
//==========================================================================//





prop = prop_pma()
prop['omega_hat'] = 0.9696
// prop['omega_hat'] = 1

var rho_eff100 = 510 // effective density
var m100 = rho_eff100 * (pi * Math.pow(100e-9, 3) / 6) // effective density @ 1 nm
prop['Dm'] = 2.48
prop['m0'] = m100 * Math.pow((1/100), prop['Dm']) // adjust mass-mobility relation parameters
var m_star = 0.01e-18

console.log('prop = ')
console.log(prop)
console.log(' ')

console.log('m_star = ')
console.log(m_star)
console.log(' ')

__left0__ = get_setpoint(prop, 'm_star', m_star, 'Rm', 5)
var sp = __left0__[0]

console.log('sp = ')
console.log(sp)
console.log(' ')

var mr_vec = linspace(1e-5, 3.7, 601)
var m_vec = mr_vec.map(function(x) {
  return x * m_star;
});

var d = m_vec.map(function(x) {
  return (Math.pow(x / prop['m0'], 1 / prop['Dm']) * 1e-9);
})

var z_vec = [1,2,3]

var Lambda_1C = parse_fun(sp, m_vec, d, prop, tfer_1C)
var Lambda_1C_diff = parse_fun(sp, m_vec, d, prop, tfer_1C_diff)
var Lambda_1S = parse_fun(sp, m_vec, d, prop, tfer_1S)
if (prop['omega_hat'] == 1) {
  var Lambda_W1 = parse_fun(sp, m_vec, d, prop, tfer_W1)
  var Lambda_W1_diff = parse_fun(sp, m_vec, d, prop, tfer_W1_diff)
}







//------------------------------------------------------------------------//
// GENERATE PLOTS

// write rho_eff100 and Dm to HTML outputs
document.getElementById('rhonum').value = rho_eff100
document.getElementById('Dmnum').value = prop['Dm']

// read resolution and mass setpoint sliders
var Rmvals = [0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15]
function displayRmval(val) {
  document.getElementById('Rmval').value = Rmvals[val - 1];
}
displayRmval(document.getElementById('RmSlider').value)

var mvals = [5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3,
  0.01, 0.02, 0.05, 0.1, 0.2, 0.5,
  1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
function displaymval(val) {
  document.getElementById('mval').value = mvals[val - 1];
}
displaymval(document.getElementById('mSlider').value)


// for legend
var margin_legend = {
  top: 0,
  right: 50,
  bottom: 0,
  left: 60
}
var svg_legend = d3.select("#my_legend")
  .append("svg")
  .attr("width", 250 + margin_legend.left + margin_legend.right)
  .attr("height", 98)
  .append("g");

// legend for lines
svg_legend.append("text")
  .attr("x", 25).attr("y", 10)
  .text("Case 1S (Ehara et al., Olfert and Collings)")
  .attr("alignment-baseline", "middle")
d1s = [{
  x: 5,
  y: 9
}, {
  x: 20,
  y: 9
}]
svg_legend.append("path")
  .datum(d1s)
  .attr("fill", "none")
  .attr("stroke", "#2525C6")
  .attr("stroke-width", 2)
  .attr("d", d3.line()
    .x(function(d) {
      return d.x;
    })
    .y(function(d) {
      return d.y;
    })
  )
svg_legend.append("text")
  .attr("x", 25).attr("y", 30)
  .text("Case 1C (Recommended over Case 1S)")
  .attr("alignment-baseline", "middle")
d1c = [{
  x: 5,
  y: 29
}, {
  x: 20,
  y: 29
}]
svg_legend.append("path")
  .datum(d1c)
  .attr("fill", "none")
  .attr("stroke", "#FFBE0B")
  .attr("stroke-width", 2)
  .attr("d", d3.line()
    .x(function(d) {
      return d.x;
    })
    .y(function(d) {
      return d.y;
    })
  )
svg_legend.append("text")
  .attr("x", 25).attr("y", 50)
  .text("Case 1C + Diffusion").attr("alignment-baseline", "middle")
d1c_diff = [{
  x: 5,
  y: 49
}, {
  x: 20,
  y: 49
}]
svg_legend.append("path")
  .datum(d1c_diff)
  .attr("fill", "none")
  .attr("stroke", "#222222")
  .attr('stroke-dasharray', "4 2")
  .attr("stroke-width", 2)
  .attr("d", d3.line()
    .x(function(d) {
      return d.x;
    })
    .y(function(d) {
      return d.y;
    })
  )
svg_legend.append("text")
  .attr("x", 25).attr("y", 70)
  .text("Case W1 (Only when ω2/ω1 = 1, where it is exact)").attr("alignment-baseline", "middle")
d1c_diff = [{
  x: 5,
  y: 69
}, {
  x: 20,
  y: 69
}]
svg_legend.append("path")
  .datum(d1c_diff)
  .attr("fill", "none")
  .attr("stroke", "#d64161")
  .attr("stroke-width", 1)
  .attr("d", d3.line()
    .x(function(d) {
      return d.x;
    })
    .y(function(d) {
      return d.y;
    })
  )
svg_legend.append("text")
  .attr("x", 25).attr("y", 90)
  .text("Case W1 + Diffusion (Only when ω2/ω1 = 1)").attr("alignment-baseline", "middle")
d1c_diff = [{
  x: 5,
  y: 89
}, {
  x: 20,
  y: 89
}]
svg_legend.append("path")
  .datum(d1c_diff)
  .attr("fill", "none")
  .attr("stroke", "#d64161")
  .attr('stroke-dasharray', "4 2")
  .attr("stroke-width", 1)
  .attr("d", d3.line()
    .x(function(d) {
      return d.x;
    })
    .y(function(d) {
      return d.y;
    })
  )






// set the dimensions and margins of the graph
var margin = {
    top: 35,
    right: 1.5,
    bottom: 50,
    left: 45
  },
  width = 750 - margin.left - margin.right,
  height = 350 - margin.top - margin.bottom;

// append the svg object to the body of the page
var svg = d3.select("#my_dataviz")
  .append("svg")
  .attr("width", width + margin.left + margin.right)
  .attr("height", height + margin.top + margin.bottom)
  .append("g")
  .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

// Add X axis
var x = d3.scaleLinear()
  .domain([0, 3.7])
  .range([0, width]);
var xAxis = svg.append("g")
  .attr("transform", "translate(0," + height + ")")
  .call(d3.axisBottom(x).ticks(5));
var xAxis2 = svg.append("g")
  .call(d3.axisTop(x).ticks(5));

// Add Y axis
var y = d3.scaleLinear()
  .domain([-0.05, 1.8])
  .range([height, 0]);
svg.append("g")
  .call(d3.axisLeft(y).ticks(5));
svg.append("g")
  .attr("transform", "translate(" + width + ",0)")
  .call(d3.axisRight(y).ticks(0))

//-- Add axis labels --//
// Add X axis label:
svg.append("text")
  .attr("text-anchor", "middle")
  .attr('x', width / 2)
  .attr('y', height + 35)
  .text("Particle mass over setpoint mass, m/m*");

// Y axis label:
svg.append("text")
  .attr("text-anchor", "middle")
  .attr('transform', 'translate(-35,' + height / 2 + ')rotate(-90)')
  .text("Transfer function, Ω")


var data = [];
for (ii = 0; ii < m_vec.length; ii++) {
  if (prop['omega_hat'] == 1) {
    data.push({
      x: mr_vec[ii],
      yc: Lambda_1C[ii],
      ys: Lambda_1S[ii],
      yc_diff: Lambda_1C_diff[ii],
      yw1: Lambda_W1[ii],
      yw1_diff: Lambda_W1_diff[ii]
    })
  } else {
    data.push({
      x: mr_vec[ii],
      yc: Lambda_1C[ii],
      ys: Lambda_1S[ii],
      yc_diff: Lambda_1C_diff[ii]
    })
  }
}

// generate plot
svg.append("path")
  .datum(data)
  .attr("id", "l1c")
  .attr("fill", "none")
  .attr("stroke", "#FFBE0B")
  .attr("stroke-width", 2)
  .attr("d", d3.line()
    .x(function(d) {
      return x(d.x)
    })
    .y(function(d) {
      return y(d.yc)
    })
  )

svg.append("path")
  .datum(data)
  .attr("id", "l1s")
  .attr("fill", "none")
  .attr("stroke", "#2525C6")
  .attr("stroke-width", 2)
  .attr("d", d3.line()
    .x(function(d) {
      return x(d.x)
    })
    .y(function(d) {
      return y(d.ys)
    })
  )

svg.append("path")
  .datum(data)
  .attr("id", "l1cd")
  .attr("fill", "none")
  .attr("stroke", "#222222")
  .attr('stroke-dasharray', "4 2")
  .attr("stroke-width", 2)
  .attr("d", d3.line()
    .x(function(d) {
      return x(d.x)
    })
    .y(function(d) {
      return y(d.yc_diff)
    })
  )

if (prop['omega_hat'] == 1) {
  svg.append("path")
    .datum(data)
    .attr("id", "lw1")
    .attr("fill", "none")
    .attr("stroke", "#d64161")
    .attr("stroke-width", 1)
    .attr("d", d3.line()
      .x(function(d) {
        return x(d.x)
      })
      .y(function(d) {
        return y(d.yw1)
      })
    )

  svg.append("path")
    .datum(data)
    .attr("id", "lw1d")
    .attr("fill", "none")
    .attr("stroke", "#d64161")
    .attr('stroke-dasharray', "4 2")
    .attr("stroke-width", 1)
    .attr("d", d3.line()
      .x(function(d) {
        return x(d.x)
      })
      .y(function(d) {
        return y(d.yw1_diff)
      })
    )
}


// adjust plot based on controls -----------------------------------------//
// control for resolution
d3.select("#RmSlider").on("change", function(d) {
  val = this.value
  Rm = Rmvals[val - 1]
  m_star = sp['m_star']
  updateData(Rm, m_star, prop)
})
// control for mass setpoint
d3.select("#mSlider").on("change", function(d) {
  val = this.value
  m_star = mvals[val - 1] / 1e18 // include conversion to kg
  Rm = sp['Rm']
  updateData(Rm, m_star, prop)
})
// control for flow rate
d3.select("#Qnum").on("change", function(d) {
  val = this.value
  prop['Q'] = val / 1000 / 60
  prop['v_bar'] = prop['Q'] / prop['A']

  Rm = sp['Rm']
  m_star = sp['m_star']
  updateData(Rm, m_star, prop)
})
// control for effective density
d3.select("#rhonum").on("change", function(d) {
  val = this.value
  rho_eff100 = val // effective density @ 100 nm read from control

  Rm = sp['Rm']
  m_star = sp['m_star']
  updateData(Rm, m_star, prop)
})
// control for mass-mobility exponent control
d3.select("#Dmnum").on("change", function(d) {
  val = this.value
  prop['Dm'] = val

  Rm = sp['Rm']
  m_star = sp['m_star']
  updateData(Rm, m_star, prop)
})
// control for omega ratio
d3.select("#omegahnum").on("change", function(d) {
  val = this.value
  prop['omega_hat'] = val

  Rm = sp['Rm']
  m_star = sp['m_star']
  updateData(Rm, m_star, prop)
})

function updateZ(data) {
  z_vec = [] // re-initialize array of integer charge states

  // For each check box:
  d3.selectAll(".cbZ").each(function(d) {
    cb = d3.select(this)
    if (cb.property("checked")) {
      z_vec.push(cb.property("value"))
    }
  })

  Rm = sp['Rm']
  m_star = sp['m_star']
  updateData(Rm, m_star, prop)
}
d3.selectAll(".cbZ").on("change", updateZ); // when button changes, run function

//------------------------------------------------------------------------//
// generic data updater
function updateData(Rm, m_star, prop) {
  m100 = rho_eff100 * (pi * Math.pow(100e-9, 3) / 6) // effective density @ 1 nm
  prop['m0'] = m100 * Math.pow((1/100), prop['Dm']) // adjust mass-mobility relation parameters

  m_vec = mr_vec.map(function(x) {
    return x * m_star;
  }) // gets points at which to evaluate the transfer function
  d = m_vec.map(function(x) {
    return (Math.pow(x / prop['m0'], 1 / prop['Dm']) * 1e-9);
  }) // gets new mobility diameters using mass-mobility relation

  // generate a new setpoint
  __left0__ = get_setpoint(prop, 'm_star', m_star, 'Rm', Rm)
  sp = __left0__[0]

  // evaulate transfer functions at new conditions
  var Lambda_1C = parse_fun(sp, m_vec, d, prop, tfer_1C)
  var Lambda_1C_diff = parse_fun(sp, m_vec, d, prop, tfer_1C_diff)
  var Lambda_1S = parse_fun(sp, m_vec, d, prop, tfer_1S)
  if (prop['omega_hat'] == 1) {
    var Lambda_W1 = parse_fun(sp, m_vec, d, prop, tfer_W1)
    var Lambda_W1_diff = parse_fun(sp, m_vec, d, prop, tfer_W1_diff)
  }

  // generate data vector to be used in updating the plot
  var data = []
  for (ii = 0; ii < m_vec.length; ii++) {
    if (prop['omega_hat'] == 1) {
      data.push({
        x: mr_vec[ii],
        yc: Lambda_1C[ii],
        ys: Lambda_1S[ii],
        yc_diff: Lambda_1C_diff[ii],
        yw1: Lambda_W1[ii],
        yw1_diff: Lambda_W1_diff[ii]
      })
    } else {
      data.push({
        x: mr_vec[ii],
        yc: Lambda_1C[ii],
        ys: Lambda_1S[ii],
        yc_diff: Lambda_1C_diff[ii]
      })
    }
  }

  // send to generiv plot updater defined below
  updatePlot(data)

  // run on update to display control values in outputs on HTML page
  document.getElementById('Vval').value = sp['V'].toPrecision(3);
  document.getElementById('Wval').value = sp['omega'].toPrecision(4);
  document.getElementById('dmval1').value =
    (Math.pow(sp['m_star'] / prop['m0'], 1 / prop['Dm'])).toPrecision(4);
  document.getElementById('dmval2').value =
    (Math.pow(2 * sp['m_star'] / prop['m0'], 1 / prop['Dm'])).toPrecision(4);
}

// run initallly to get control values and display in outputs on HTML page
document.getElementById('Vval').value = sp['V'].toPrecision(3);
document.getElementById('Wval').value = sp['omega'].toPrecision(4);
document.getElementById('dmval1').value =
  (Math.pow(sp['m_star'] / prop['m0'], 1 / prop['Dm'])).toPrecision(4);
document.getElementById('dmval2').value =
  (Math.pow(2 * sp['m_star'] / prop['m0'], 1 / prop['Dm'])).toPrecision(4);
document.getElementById('omegahnum').value = prop['omega_hat'];
//------------------------------------------------------------------------//

// a generic function that updates the chart -----------------------------//
function updatePlot(data) {
  // consider 1C case
  d3.select("#l1c")
    .datum(data)
    .transition()
    .attr("d", d3.line()
      .x(function(d) {
        return x(d.x)
      })
      .y(function(d) {
        return y(d.yc)
      })
    )

  // consider 1S case
  d3.select("#l1s")
    .datum(data)
    .transition()
    .attr("d", d3.line()
      .x(function(d) {
        return x(d.x)
      })
      .y(function(d) {
        return y(d.ys)
      })
    )

  // consider 1C diffusing case
  d3.select("#l1cd")
    .datum(data)
    .transition()
    .attr("d", d3.line()
      .x(function(d) {
        return x(d.x)
      })
      .y(function(d) {
        return y(d.yc_diff)
      })
    )

  // consider the W1 case
  // including removing line when ω2/ω1 is not unity or
  // adding line when ω2/ω1 is changed to unity
  if (prop['omega_hat'] == 1) { // add line if w2/w1 is changed to unity
    if (d3.select("#lw1").empty()) {
      svg.append("path")
        .datum(data)
        .attr("id", "lw1")
        .attr("fill", "none")
        .attr("stroke", "#d64161")
        .attr("stroke-width", 1)
        .attr("d", d3.line()
          .x(function(d) {
            return x(d.x)
          })
          .y(function(d) {
            return y(d.yw1)
          })
        )
      svg.append("path")
        .datum(data)
        .attr("id", "lw1d")
        .attr("fill", "none")
        .attr("stroke", "#d64161")
        .attr('stroke-dasharray', "4 2")
        .attr("stroke-width", 1)
        .attr("d", d3.line()
          .x(function(d) {
            return x(d.x)
          })
          .y(function(d) {
            return y(d.yw1_diff)
          })
        )
    } else { // adjust line if w2/w1 remains unity
      d3.select("#lw1")
        .datum(data)
        .transition()
        .attr("d", d3.line()
          .x(function(d) {
            return x(d.x)
          })
          .y(function(d) {
            return y(d.yw1)
          })
        )
      d3.select("#lw1d")
        .datum(data)
        .transition()
        .attr("d", d3.line()
          .x(function(d) {
            return x(d.x)
          })
          .y(function(d) {
            return y(d.yw1_diff)
          })
        )
    }
  } else { // remove line if w2/w1 is no longer unity
    if (!(d3.select("#lw1").empty())) {
      d3.select("#lw1").remove();
      d3.select("#lw1d").remove();
    }
  }
}


// add labels for charge states
svg.append("circle")
  .attr("cx", function(d) { return x(1) })
  .attr("cy", -18)
  .attr("r", 10)
  .style("fill", "#fff2cc")
  .attr("stroke", "black")
  .attr("stroke-width", 1.2)
svg.append("text")
  .attr("text-anchor", "middle")
  .attr('transform', 'translate(' + x(1) + ',-14)')
  .style("font-size", "11px")
  .text("+1")
svg.append("circle")
  .attr("cx", function(d) { return x(2) })
  .attr("cy", -18)
  .attr("r", 10)
  .style("fill", "#fff2cc")
  .attr("stroke", "black")
  .attr("stroke-width", 1.2)
svg.append("text")
  .attr("text-anchor", "middle")
  .attr('transform', 'translate(' + x(2) + ',-14)')
  .style("font-size", "11px")
  .text("+2")
svg.append("circle")
  .attr("cx", function(d) { return x(3) })
  .attr("cy", -18)
  .attr("r", 10)
  .style("fill", "#fff2cc")
  .attr("stroke", "black")
  .attr("stroke-width", 1.2)
svg.append("text")
  .attr("text-anchor", "middle")
  .attr('transform', 'translate(' + x(3) + ',-14)')
  .style("font-size", "11px")
  .text("+3")
svg.append("text")
  .attr("text-anchor", "left")
  .attr('transform', 'translate(' + 40 + ',-13)')
  .style("font-size", "13px")
  .text("Integer charge state > ")

//------------------------------------------------------------------------//
//END PLOT ---------------------------------------------------------------//
