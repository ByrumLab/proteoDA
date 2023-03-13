function makeVegaTooltip(columns)
{
    // generate tooltip object for embedding in spec
    let tooltipString = "{";
    // ignore the internal_id_for_brushing data, otherwise all
    // data passed from R goes into the tooltip
    columns.forEach( x => {if (x !== "internal_id_for_brushing") {tooltipString += `'${x}':datum['${x}'],`}});
    tooltipString += "}";
    var tooltip = { "signal" : tooltipString };
    return tooltip;
}
