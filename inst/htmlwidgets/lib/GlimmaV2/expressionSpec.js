function createExpressionSpec(width, height, expColumns, sampleColours, samples, numUniqueGroups)
{

    /* need an empty signal if sample.cols argument has been supplied */
    let samplecols_signal = { "name": "samplecols_active" };

    /* must match counts term in processExpression */
    expColumns.push("normalized intensity");
    let tooltip = makeVegaTooltip(expColumns);
    let tooltip_mean = makeVegaTooltip(["group", "mean"]);

    return {
        "$schema": "https://vega.github.io/schema/vega/v5.json",
        "width": width*0.40,
        "height": height*0.35,
        "padding": {"left": 0, "top": 0, "right": 0, "bottom": 10},
        "autosize": {"type": "fit", "resize": true},
        "title": { "text": {"signal": "title_signal" }},
        "signals":
                [
                    {
                        "name": "title_signal",
                        "value": ""
                    },
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
                        "update": " (ylim_min === null || ylim_min == \"\") ? null : (ylim_min > min_count) ? null : ylim_min"
                    },
                    {
                        "name": "max_y",
                        "update": " (ylim_max < max_count) ? null : ylim_max"
                    },
                    samplecols_signal,
                    {
                      "name": "mean_linewidth",
                      "value":  ((width*0.4)/(numUniqueGroups + 1)*0.6)*((width*0.4)/(numUniqueGroups + 1)*0.6)
                    }
                ],
        "data": [
          {"name": "table"},
          {"name": "means",
           "source" : "table",
           "transform" :
             [{
              "type": "aggregate",
              "groupby": ["group"],
              "fields": ["normalized intensity"],
              "ops": ["mean"],
              "as": ["mean"]
             }]
          }
        ],
        "scales":
        [
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
        "axes":
        [
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
          {
            "name": "med_lines",
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
