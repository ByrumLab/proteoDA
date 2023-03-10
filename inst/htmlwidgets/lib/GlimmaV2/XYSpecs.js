// Returns the vega spec for the left side (XY) plot
function createXYSpec(xyData, xyTable, width, height)
{
  // Make the tooltip string first using the JS function
  // in lib/GlimmaV2/makeTooltip.js
  var tooltip = makeVegaTooltip(xyData.cols);

  return {
    "$schema": "https://vega.github.io/schema/vega/v5.json",
    "description": "Testing ground for GlimmaV2",
    "width": width * 0.5,
    "height": height * 0.35,
    "padding": {"left": 0, "top": 0, "right": 0, "bottom": 10},
    "autosize": {"type": "fit", "resize": true},
    "signals":
      [
        {
          "name": "click", "value": null,
          "on": [ {"events": "mousedown", "update": "[datum, now()]" } ]
        },
        // Some custom signals
        // The two main ones are the plot_type and pval_type,
        // which are bound to selection boxes using vega's selection binding
        // Once plot type and pval type are chosen, those signals are passed
        // down to other signals which output the appropriate columns to use for
        // the axes, axes titles, etc.
        {
          "name": "plot_type",
          "value": "volcano",
          "bind": {
            "input": "select",
            "options": ["volcano", "MD"],
            "name": "Plot type: ",
            "style": "width: 200px; margin-bottom: 10px;"
          }
        },
        {
          "name": "pval_type",
          "value": "adjusted",
          "bind": {
            "input": "select",
            "options": ["raw", "adjusted"],
            "name": "P-value type: ",
            "style": "width: 200px; margin-bottom: 10px;"
          }
        },
        // Set x axis data and title based on plot type
        // Only two options
        {
          "name": "x_axis",
          "value": "logFC",
          "on": [
            {
             "events": [
                {"signal": "pval_type"},
                {"signal": "plot_type"}
                ],
              "update": "(plot_type == \"MD\" ? \"average_intensity\" : \"logFC\")"
            }
          ]
        },
        {
          "name": "x_axis_title",
          "value": "log2 fold-change (logFC)",
          "on": [
            {
             "events": [
                {"signal": "pval_type"},
                {"signal": "plot_type"}
                ],
              "update": "(plot_type == \"MD\" ? \"average intensity\" : \"log2 fold-change (logFC)\")"
            }
          ]
        },
        // Setting y axis is slightly more complicated, with three options.
        // vega has slightly limited syntax, can't use a complicated function,
        // just use nested ternary opertors.
        {
          "name": "y_axis",
          "value": "negLog10adjP",
           "on": [
            {
              "events": [
                {"signal": "pval_type"},
                {"signal": "plot_type"}
                ],
              "update": "(plot_type == \"MD\" ? \"logFC\" : (pval_type == \"raw\" ? \"negLog10rawP\" : \"negLog10adjP\"))"
            }
            ]
        },
                {
          "name": "y_axis_title",
          "value": "-log10(adjusted P)",
           "on": [
            {
              "events": [
                {"signal": "pval_type"},
                {"signal": "plot_type"}
                ],
              "update": "(plot_type == \"MD\" ? \"log2 fold-change (logFC)\" : (pval_type == \"raw\" ? \"-log10(raw P)\" : \"-log10(adjusted P)\"))"
            }
            ]
        },
        // point color scale just based on p-value type
        {
          "name": "point_color",
          "value": "sig.FDR.fct",
           "on": [
            {
              "events": [
                {"signal": "pval_type"},
                {"signal": "plot_type"}
                ],
              "update": "(pval_type == \"raw\" ? \"sig.pval.fct\" : \"sig.FDR.fct\")"
            }
            ]
        }
      ],
    "data":
      [
        {
          "name": "source",
          "values": xyTable,
          "transform": [{
            "type": "formula",
            "expr": "datum.x",
            "as": "tooltip"
          }]
        },
        { "name": "selected_points" }
      ],
    "scales": [
      {
        "name": "x",
        "type": "linear",
        "round": true,
        "nice": true,
        "zero": false,
        "domain": { "data": "source", "field": {"signal" : "x_axis"}},
        "range": "width"
      },
      {
        "name": "y",
        "type": "linear",
        "round": true,
        "nice": true,
        "zero": false,
        "domain": {
          "data": "source",
          "field": {"signal" : "y_axis"}
        },
        "range": "height"
      },
      {
        "name": "colour_scale",
        "type": "ordinal",
        "domain": ["downReg", "nonDE", "upReg"],
        "range": xyData.statusColours
      }
    ],
    "legends": [
      {
        "fill": "colour_scale",
        "title": "Status",
        "symbolStrokeColor": "black",
        "symbolStrokeWidth": 1,
        "symbolOpacity": 0.7,
        "symbolType": "circle"
      }
    ],
    "axes" : [
      {
        "scale": "x",
        "grid": true,
        "domain": false,
        "orient": "bottom",
        "tickCount": 5,
        "title": {"signal" : "x_axis_title"}
      },
      {
        "scale": "y",
        "grid": true,
        "domain": false,
        "orient": "left",
        "titlePadding": 5,
        "title": {"signal" : "y_axis_title"}
      }
    ],
    "marks": [
      {
        "name": "marks",
        "type": "symbol",
        "from": { "data": "source" },
        "encode": {
          "update": {
            "x": { "scale": "x", "field": {"signal" : "x_axis"} },
            "y": {
              "scale": "y",
              "field": {"signal" : "y_axis"}
            },
            "shape": "circle",
            "size" : {"value": 25},
            "opacity": {"value": 0.65},
            "fill": { "scale": "colour_scale", "field": {"signal": "point_color"} },
            "strokeWidth": {"value": 1},
            "stroke": {"value": "transparent"},
            "tooltip": tooltip
          }
        }
      },
      // overlaying selected points
      // with a larger mark when selected
      {
        "name": "selected_marks",
        "type": "symbol",
        "from": { "data": "selected_points" },
        "encode": {
          "update": {
            "x": { "scale": "x", "field": {"signal" : "x_axis"}},
            "y": {
              "scale": "y",
              "field": {"signal" : "y_axis"}
            },
            "shape": "circle",
            "size": {"value": 120},
            "fill": { "scale": "colour_scale", "field":  {"signal": "point_color"}  },
            "strokeWidth": { "value": 1 },
            "stroke": { "value": "black" },
            "opacity": { "value": 1 },
            "tooltip": tooltip
          }
        }
      }
    ]
  };
}
