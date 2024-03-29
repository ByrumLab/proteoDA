// Returns the vega spec for the right side (abundance) plot
// called "expression" or "exp" across most of the glimma code,
// which was originally for RNAseq stuff
function createExpressionSpec(width, height, expColumns, sampleColours, numUniqueGroups)
{
  // Add the normalized intensity data to the expColumns,
  // so it gets included in the tooltip
  expColumns.push("normalized intensity");
  // Make the tooltip strings first using the JS function
  // in lib/GlimmaV2/makeTooltip.js
  let tooltip = makeVegaTooltip(expColumns);
  let tooltip_mean = makeVegaTooltip(["group", "mean"]);

  return {
    "$schema": "https://vega.github.io/schema/vega/v5.json",
    "width": width*0.40,
    "height": height*0.35,
    "padding": {"left": 0, "top": 0, "right": 0, "bottom": 10},
    "autosize": {"type": "fit", "resize": true},
    "title": { "text": {"signal": "title_signal" }},
    "signals": [
      {
        "name": "title_signal",
        "value": ""
      },
      // Users can set custom ylims bound to vega
      // number selection fields
      {
        "name": "ylim_min",
        "value": null,
        "bind": {
          "input": "number",
          "class": "ylim_min"
        }
      },
      {
        "name": "ylim_max",
        "value": null,
        "bind": {
          "input": "number",
          "class": "ylim_max"
        }
      },
      // These signals then check the min and max
      // values of the data coming in and pass that on to
      // the domain for the y axis
      {
        "name": "min_count",
        "value": null
      },
      {
        "name": "max_count",
        "value": 0
      },
      {
        "name": "min_y",
        "update": "(ylim_min === null || ylim_min == \"\") ? null : (ylim_min > min_count) ? null : ylim_min"
      },
      {
        "name": "max_y",
        "update": "(ylim_max < max_count) ? null : ylim_max"
      },
      // length/width of the mean line is the square of the desired length in pixels.
      // calculation is similar to the calcs done for offset in
      // glimmaXY.js, using 60% of the pixels alloted to each group
      {
        "name": "mean_linewidth",
        "value":  ((width*0.4)/(numUniqueGroups + 1)*0.6)*((width*0.4)/(numUniqueGroups + 1)*0.6)
      }
    ],
    "data": [
      // main data comes from the table signal, passed in from the JS code
      // in glimmaXY.js
      {"name": "table"},
      // Then, from the incomign data,
      // calculate the means using vega transfroms and aggregates
      {"name": "means",
       "source" : "table",
       "transform" : [
         {
           "type": "aggregate",
           "groupby": ["group"],
           "fields": ["normalized intensity"],
           "ops": ["mean"],
           "as": ["mean"]
         }
       ]
      }
    ],
    "scales": [
      {
        "name": "x",
        "type": "point",
        "padding": 1,
        "domain": {"data": "table", "field": "group"},
        "range": "width"
      },
      {
        "name": "y",
        "domain": {"data": "table", "field": "normalized intensity"},
        "range": "height",
        "zero": false,
        "domainMin": {"signal": "min_y"},
        "domainMax": {"signal": "max_y"}
      },
      {
        "name": "color",
        "type": "ordinal",
        "domain": { "data": "table", "field": "group" },
        "range": sampleColours
      }
    ],
    "axes": [
      {
        "scale": "x",
        "orient": "bottom",
        "title": "group",
        "labelAngle": -45,
        "labelAlign": "right",
        "labelOffset": -3
      },
      {
        "scale": "y",
        "grid": true,
        "orient": "left",
        "titlePadding": 5,
        "title": "normalized intensity"
      }
    ],
    "marks": [
      // The main marks of intensity
      {
        "name": "marks",
        "type": "symbol",
        "from": {"data": "table"},
        "encode": {
          "update": {
            "x": {"scale": "x", "field": "group", "offset": {"field": "offset"}},
            "y": {"scale": "y", "field": "normalized intensity"},
            "shape": {"value": "circle"},
            "fill": { "scale": "color", "field": "group"},
            "strokeWidth": {"value": 1},
            "opacity": {"value": 0.8},
            "size": {"value": 100},
            "tooltip": tooltip
          }
        }
      },
      // The mean lines
      {
        "name": "mean_lines",
        "type": "symbol",
        "from": {"data": "means"},
        "encode": {
          "update": {
            "x": {"scale": "x", "field": "group"},
            "y": {"scale": "y", "field": "mean"},
            "shape": {"value": "stroke"},
            "stroke": {"scale" : "color", "field": "group"},
            "strokeWidth": {"value": 3.5},
            "opacity": {"value": 1},
            "size": {"signal": "mean_linewidth"},
            "tooltip": tooltip_mean
          }
        }
      }
    ]
  };
}
