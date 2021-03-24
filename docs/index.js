prop = prop_pma()
prop['omega_hat'] = 0.9696
// prop['omega_hat'] = 1

var rho_eff100 = 510 // effective density
var m100 = rho_eff100 * (pi * Math.pow(100e-9, 3) / 6) // effective density @ 1 nm
prop['Dm'] = 2.48
prop['m0'] = m100 * Math.pow((1 / 100), prop['Dm']) // adjust mass-mobility relation parameters
var m_star = 0.01e-18

console.log('prop = ')
console.log(prop)
console.log(' ')

console.log('m_star = ')
console.log(m_star)
console.log(' ')

sp_var1 = 'm_star';
sp_var2 = 'Rm';
__left0__ = get_setpoint(prop, sp_var1, m_star, sp_var2, 5)
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

var z_vec = [1, 2, 3]

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

var mvals = [5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3,
  0.01, 0.02, 0.05, 0.1, 0.2, 0.5,
  1, 2, 5, 10, 20, 50, 100, 200, 500, 1000
]

// Define color scheme.
var colors = ["#2525C6", "#FFBE0B", "#222222", "#D64161"]; // original
var colors = ["#03C9A5", "#1E3778", "#222222", "#DD393A"]; // new
var colors = ["#1E2978", "#48DEB1", "#222222", "#DD393A"]; // new

function displaymval() {
  document.getElementById('mval').innerHTML = mvals[document.getElementById('mSlider').value - 1];
}
displaymval()



var $container = $('#my_dataviz'),
  width_a = 0.95 * Math.min($container.width(), 870),
  height_a = $container.height()


// for legend
var margin_legend = {
  top: 0,
  right: 50,
  bottom: 0,
  left: 60
}
var svg_legend = d3.select("#my_legend")
  .append("svg")
  .attr("width", width_a - margin_legend.left - margin_legend.right)
  .attr("height", 110)
  .append("g");

// legend for lines
svg_legend.append("text")
  .attr("x", 25).attr("y", 12)
  .text("Case 1S (Ehara et al., Olfert and Collings)")
  .attr("alignment-baseline", "middle");

d1s = [{
  x: 0,
  y: 11
}, {
  x: 18,
  y: 11
}]
svg_legend.append("path")
  .datum(d1s)
  .attr("fill", "none")
  .attr("stroke", colors[0])
  .attr("stroke-width", 2.75)
  .attr("d", d3.line()
    .x(function(d) {
      return d.x;
    })
    .y(function(d) {
      return d.y;
    })
  )
svg_legend.append("text")
  .attr("x", 25).attr("y", 34)
  .text("Case 1C (Recommended over Case 1S)")
  .attr("alignment-baseline", "middle")
d1c = [{
  x: 0,
  y: 33
}, {
  x: 18,
  y: 33
}]
svg_legend.append("path")
  .datum(d1c)
  .attr("fill", "none")
  .attr("stroke", colors[1])
  .attr("stroke-width", 2.75)
  .attr("d", d3.line()
    .x(function(d) {
      return d.x;
    })
    .y(function(d) {
      return d.y;
    })
  )
svg_legend.append("text")
  .attr("x", 25).attr("y", 56)
  .text("Case 1C + Diffusion").attr("alignment-baseline", "middle")
d1c_diff = [{
  x: 0,
  y: 54
}, {
  x: 18,
  y: 54
}]
svg_legend.append("path")
  .datum(d1c_diff)
  .attr("fill", "none")
  .attr("stroke", colors[2])
  .attr('stroke-dasharray', "4 2")
  .attr("stroke-width", 2.75)
  .attr("d", d3.line()
    .x(function(d) {
      return d.x;
    })
    .y(function(d) {
      return d.y;
    })
  )
svg_legend.append("text")
  .attr("x", 25).attr("y", 78)
  .text("Case W1 (Only when ω2/ω1 = 1, where it is exact)").attr("alignment-baseline", "middle")
d1c_diff = [{
  x: 0,
  y: 77
}, {
  x: 18,
  y: 77
}]
svg_legend.append("path")
  .datum(d1c_diff)
  .attr("fill", "none")
  .attr("stroke", colors[3])
  .attr("stroke-width", 1.0)
  .attr("d", d3.line()
    .x(function(d) {
      return d.x;
    })
    .y(function(d) {
      return d.y;
    })
  )
svg_legend.append("text")
  .attr("x", 25).attr("y", 100)
  .text("Case W1 + Diffusion (Only when ω2/ω1 = 1)").attr("alignment-baseline", "middle")
d1c_diff = [{
  x: 0,
  y: 99
}, {
  x: 18,
  y: 99
}]
svg_legend.append("path")
  .datum(d1c_diff)
  .attr("fill", "none")
  .attr("stroke", colors[3])
  .attr('stroke-dasharray', "4 2")
  .attr("stroke-width", 1.0)
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
    top: 0,
    right: 1.5,
    bottom: 50,
    left: 65
  },
  width = width_a - margin.left - margin.right,
  height = 350 - margin.top - margin.bottom;

// append the svg object to the body of the page
var svg = d3.select("#my_dataviz")
  .append("svg")
  .attr("width", width + margin.left + margin.right)
  .attr("height", height + margin.top + margin.bottom)
  .append("g")
  .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

//-- Add background rectangle --//
svg.append("rect")
  .attr("width", "100%")
  .attr("height", height)
  .attr("fill", "#FFF");

// Add X axis
var x = d3.scaleLinear()
  .domain([0, 3.7])
  .range([0, width]);
var xAxis = svg.append("g")
  .attr("transform", "translate(0," + height + ")")
  .attr("class", "axis")
  .call(d3.axisBottom(x).ticks(5));
var xAxis2 = svg.append("g")
  .attr("class", "axis")
  .call(d3.axisTop(x).ticks(0));

// Add Y axis
var yMax = 1.6,
  yMin = -0.05;

var y = d3.scaleLinear()
  .domain([yMin, yMax])
  .range([height, 0]);
var yAxis = svg.append("g")
  .attr("class", "axis")
  .call(d3.axisLeft(y).ticks(5));
var yAxis2 = svg.append("g")
  .attr("transform", "translate(" + width + ",0)")
  .attr("class", "axis")
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
  .attr('transform', 'translate(-40,' + height / 2 + ')rotate(-90)')
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
  .attr("stroke", colors[1])
  .attr("stroke-width", 2.75)
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
  .attr("id", "l1cd")
  .attr("fill", "none")
  .attr("stroke", colors[2])
  .attr('stroke-dasharray', "4 2")
  .attr("stroke-width", 2.75)
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
    .attr("stroke", colors[3])
    .attr("stroke-width", 1.0)
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
    .attr("stroke", colors[3])
    .attr('stroke-dasharray', "4 2")
    .attr("stroke-width", 1.0)
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
d3.select("#RmSlider").on("change", function() {
  val1 = parseFloat(this.value)
  if (val1 <= 0) {
    console.log(val1)
    $('#RmSlider').css("border", "2px solid #EE696A");
    return;
  }
  $('#RmSlider').hover(function() {
    $(this).css("border", "2px solid #9BD6FF");
  }, function() {
    $(this).css("border", "2px solid #E3E3E3");
  });
  Rm = val1;
  m_star = sp[sp_var1]
  updateData(Rm, m_star, prop)
})
// control for mass setpoint
d3.select("#mSlider").on("change", function() {
  val = this.value
  m_star = mvals[val - 1] / 1e18 // include conversion to kg
  Rm = sp[sp_var2]
  updateData(Rm, m_star, prop)
})
// control for flow rate
d3.select("#Qnum").on("change", function() {
  val = this.value
  prop['Q'] = val / 1000 / 60
  prop['v_bar'] = prop['Q'] / prop['A']

  Rm = sp[sp_var2]
  m_star = sp[sp_var1]
  updateData(Rm, m_star, prop)
})
// control for effective density
d3.select("#rhonum").on("change", function() {
  val = this.value
  rho_eff100 = val // effective density @ 100 nm read from control

  Rm = sp[sp_var2]
  m_star = sp[sp_var1]
  updateData(Rm, m_star, prop)
})
// control for mass-mobility exponent control
d3.select("#Dmnum").on("change", function() {
  val = this.value
  prop['Dm'] = val

  Rm = sp[sp_var2]
  m_star = sp[sp_var1]
  updateData(Rm, m_star, prop)
})
// control for omega ratio
d3.select("#omegahnum").on("change", function() {
  val = this.value
  prop['omega_hat'] = val

  Rm = sp[sp_var2]
  m_star = sp[sp_var1]
  updateData(Rm, m_star, prop)
})

// control for length of classifier
d3.select("#Lnum").on("change", function() {
  val = this.value
  prop['L'] = val

  Rm = sp[sp_var2]
  m_star = sp[sp_var1]
  updateData(Rm, m_star, prop)
})

// control for r1
function afterRadiusUpdate(prop) {
  prop['del'] = (prop['r2'] - prop['r1']) / 2
  prop['rc'] = (prop['r2'] + prop['r1']) / 2
  prop['r_hat'] = prop['r1'] / prop['r2']
  prop['A'] = pi * (Math.pow(prop['r2'], 2) - Math.pow(prop['r1'], 2));
  prop['v_bar'] = prop['Q'] / prop['A'];
  return prop;
}
d3.select("#r1num").on("change", function() {
  val = this.value
  prop['r1'] = val / 100
  prop = afterRadiusUpdate(prop)
  document.getElementById('r2num').min = prop['r1'] * 100 + 0.005 // prevent overlapping radii

  Rm = sp[sp_var2]
  m_star = sp[sp_var1]
  updateData(Rm, m_star, prop)
})
// control for r2
d3.select("#r2num").on("change", function() {
  val = this.value
  prop['r2'] = val / 100
  prop = afterRadiusUpdate(prop)
  document.getElementById('r1num').max = prop['r2'] * 100 - 0.005 // prevent overlapping radii

  Rm = sp[sp_var2]
  m_star = sp[sp_var1]
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

  Rm = sp[sp_var2]
  m_star = sp[sp_var1]
  updateData(Rm, m_star, prop)
}
d3.selectAll(".cbZ").on("change", updateZ); // when button changes, run function

// for links in text to set specific conditions (e.g., jump to APM)
d3.select("#setapm").on("click", function() {
  prop['omega_hat'] = 1
  document.getElementById('omegahnum').value = 1;

  prop['Q'] = 0.5 / 1000 / 60;
  document.getElementById('Qnum').value = 0.5;

  prop['r1'] = 0.1;
  prop['r2'] = 0.103;
  prop = afterRadiusUpdate(prop)

  document.getElementById('r1num').value = 10;
  document.getElementById('r2num').value = 10.3;
  document.getElementById('r1num').max = 10.3 - 0.005;
  document.getElementById('r2num').min = 10 + 0.005;

  Rm = sp[sp_var2]
  m_star = sp[sp_var1]
  updateData(Rm, m_star, prop)
})

d3.select("#setcpma").on("click", function() {
  prop['omega_hat'] = 0.9696
  document.getElementById('omegahnum').value = 1

  prop['Q'] = 3 / 1000 / 60;
  document.getElementById('Qnum').value = 3

  prop['r1'] = 0.06;
  prop['r2'] = 0.061;
  prop = afterRadiusUpdate(prop)
  document.getElementById('r1num').value = 0.06 * 100
  document.getElementById('r2num').value = 0.061 * 100
  document.getElementById('r1num').max = 0.061 * 100 - 0.005
  document.getElementById('r2num').min = 0.06 * 100 + 0.005

  Rm = sp[sp_var2]
  m_star = sp[sp_var1]
  updateData(Rm, m_star, prop);
})

d3.select("#fCharge").on("click", function() {
  fCharge = 1 - fCharge;
  Rm = sp[sp_var2]
  m_star = sp[sp_var1]
  updateData(Rm, m_star, prop)
})

// Define variable for broader scope.
// Flag for wheter to plot Case 1S.
var f1S = false
d3.select("#f1S").on("click", function() {
  f1S = document.getElementById('f1S').checked
  Rm = sp[sp_var2]
  m_star = sp[sp_var1]
  updateData(Rm, m_star, prop)
})




//------------------------------------------------------------------------//
// Different setpoint modes.
d3.select("#sp-mode").on("change", function() {
  mode = this.value;

  switch (mode) {
    case "Mass + Resolution":
      sp_var1 = "m_star"
      sp_var2 = "Rm"
      document.getElementById('var1-name').innerHTML = 'Mass setpoint';
      document.getElementById('var2-name').innerHTML = 'Resolution';
      document.getElementById('var2-units').innerHTML = '';
      document.getElementById('RmSlider').value = sp['Rm'].toPrecision(3);
      break;

    case "Mass + Voltage":
      sp_var1 = "m_star"
      sp_var2 = "V"
      document.getElementById('var1-name').innerHTML = 'Mass setpoint';
      document.getElementById('var2-name').innerHTML = 'Voltage';
      document.getElementById('var2-units').innerHTML = 'V';
      document.getElementById('RmSlider').value = sp['V'].toPrecision(3);
      break;

    case "Mass + Angular speed":
      sp_var1 = "m_star"
      sp_var2 = "omega"
      document.getElementById('var1-name').innerHTML = 'Mass setpoint';
      document.getElementById('var2-name').innerHTML = 'Angular speed';
      document.getElementById('var2-units').innerHTML = 'rad/s';
      document.getElementById('RmSlider').value = sp['omega'].toFixed(0);
      break;

  }

  updateData(sp[sp_var2], sp[sp_var1], prop)
})



//------------------------------------------------------------------------//
// Generic data updater.
function updateData(Rm, m_star, prop) {
  console.log("var1 = " + sp_var1 + ", " + m_star)
  console.log("var2 = " + sp_var2 + ", " + Rm)
  console.log(" ")

  m100 = rho_eff100 * (pi * Math.pow(100e-9, 3) / 6) // effective density @ 1 nm
  prop['m0'] = m100 * Math.pow((1 / 100), prop['Dm']) // adjust mass-mobility relation parameters

  m_vec = mr_vec.map(function(x) {
    return x * m_star;
  }) // gets points at which to evaluate the transfer function
  d = m_vec.map(function(x) {
    return (Math.pow(x / prop['m0'], 1 / prop['Dm']) * 1e-9);
  }) // gets new mobility diameters using mass-mobility relation

  // generate a new setpoint
  __left0__ = get_setpoint(prop, sp_var1, m_star, sp_var2, Rm)
  sp = __left0__[0]

  console.log("sp = ")
  console.log(sp)
  console.log(" ")

  // evaulate transfer functions at new conditions
  var Lambda_1C = parse_fun(sp, m_vec, d, prop, tfer_1C)
  var Lambda_1C_diff = parse_fun(sp, m_vec, d, prop, tfer_1C_diff)
  var Lambda_1S = parse_fun(sp, m_vec, d, prop, tfer_1S)
  if (prop['omega_hat'] == 1) {
    var Lambda_W1 = parse_fun(sp, m_vec, d, prop, tfer_W1)
    var Lambda_W1_diff = parse_fun(sp, m_vec, d, prop, tfer_W1_diff)
  }
  yMax = 1.6 * Math.max.apply(this, Lambda_1C);
  yMin = -0.05 * Math.max.apply(this, Lambda_1C);

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

  // run on update to display control values in outputs on HTML page
  updateVals()

  // send to generiv plot updater defined below
  updatePlot(data)
}

// run initallly to get control values and display in outputs on HTML page
console.log(document.getElementById('Vval').innerHTML)
var updateVals = function() {
  document.getElementById('Vval').innerHTML = sp['V'].toPrecision(3);
  document.getElementById('Wval').innerHTML = sp['omega'].toPrecision(4);
  document.getElementById('Rmval2').innerHTML = sp['Rm'].toPrecision(3);
  document.getElementById('dmval1').innerHTML =
    (Math.pow(sp['m_star'] / prop['m0'], 1 / prop['Dm'])).toPrecision(4);
  document.getElementById('dmval2').innerHTML =
    (Math.pow(2 * sp['m_star'] / prop['m0'], 1 / prop['Dm'])).toPrecision(4);
  document.getElementById('dmval3').innerHTML =
    (Math.pow(3 * sp['m_star'] / prop['m0'], 1 / prop['Dm'])).toPrecision(4);
  document.getElementById('omegahnum').value = prop['omega_hat'];

  document.getElementById('r1num').value = prop['r1'] * 100;
  document.getElementById('r1num').max = prop['r2'] * 100 - 0.005;
  document.getElementById('r2num').value = prop['r2'] * 100;
  document.getElementById('r2num').min = prop['r1'] * 100 + 0.005;
}
updateVals()
//------------------------------------------------------------------------//

// a generic function that updates the chart -----------------------------//
function updatePlot(data) {

  y = d3.scaleLinear()
    .domain([yMin, yMax])
    .range([height, 0]);
  yAxis.call(d3.axisLeft(y).ticks(5));
  yAxis2.attr("transform", "translate(" + width + ",0)")
    .call(d3.axisRight(y).ticks(0))

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
  if (f1S) {
    if (d3.select("#l1s").empty()) {
      svg.append("path")
        .datum(data)
        .attr("id", "l1s")
        .attr("fill", "none")
        .attr("stroke", colors[0])
        .attr("stroke-width", 2.75)
        .attr("d", d3.line()
          .x(function(d) {
            return x(d.x)
          })
          .y(function(d) {
            return y(0 * d.ys)
          })
        )
    }
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
  } else { // remove line if w2/w1 is no longer unity
    if (!(d3.select("#l1s").empty())) {
      d3.select("#l1s")
        .transition()
        .attr("stroke-width", 0);

      setTimeout( function () {  d3.select("#l1s").remove(); }, 200);
    }
  }

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
        .attr("stroke", colors[3])
        .attr("stroke-width", 1.0)
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
        .attr("stroke", colors[3])
        .attr('stroke-dasharray', "4 2")
        .attr("stroke-width", 1.0)
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
updatePlot(data); // run the plot update initially



// add labels for charge states
// set the dimensions and margins of the graph
var marginl = {
    top: 30,
    right: margin.right,
    bottom: 0,
    left: margin.left
  },
  widthl = width_a - marginl.left - marginl.right,
  heightl = 32 - marginl.top - marginl.bottom;

// append the svg object to the body of the page
var svg3 = d3.select("#myz")
  .append("svg")
  .attr("width", widthl + marginl.left + marginl.right)
  .attr("height", heightl + marginl.top + marginl.bottom)
  .append("g")
  .attr("transform", "translate(" + marginl.left + "," + marginl.top + ")");
var xAxis3 = svg3.append("g")
  .call(d3.axisTop(x).tickValues([1, 2, 3]))

svg3.append("circle")
  .attr("cx", function() {
    return x(1)
  })
  .attr("cy", -18)
  .attr("r", 11)
  .style("fill", "#ededed")
  .attr("stroke", "black")
  .attr("stroke-width", 1.0)
svg3.append("text")
  .attr("text-anchor", "middle")
  .attr('transform', 'translate(' + x(1) + ',-14)')
  .style("font-size", "11px")
  .text("+1")
svg3.append("circle")
  .attr("cx", function() {
    return x(2)
  })
  .attr("cy", -18)
  .attr("r", 11)
  .style("fill", "#ededed")
  .attr("stroke", "black")
  .attr("stroke-width", 1.0)
svg3.append("text")
  .attr("text-anchor", "middle")
  .attr('transform', 'translate(' + x(2) + ',-14)')
  .style("font-size", "11px")
  .text("+2")
svg3.append("circle")
  .attr("cx", function() {
    return x(3)
  })
  .attr("cy", -18)
  .attr("r", 11)
  .style("fill", "#ededed")
  .attr("stroke", "black")
  .attr("stroke-width", 1.0)
svg3.append("text")
  .attr("text-anchor", "middle")
  .attr('transform', 'translate(' + x(3) + ',-14)')
  .style("font-size", "11px")
  .text("+3")
svg3.append("text")
  .attr("text-anchor", "left")
  .attr('transform', 'translate(' + 30 + ',-13)')
  .style("font-size", "13px")
  .text("z = ")

//------------------------------------------------------------------------//
//END PLOT ---------------------------------------------------------------//
