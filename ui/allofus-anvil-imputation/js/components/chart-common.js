/**
 * Shared styling and helpers for the validation charts (js/components/*-chart.js).
 *
 * Chart types register themselves in VALIDATION_CHART_TYPES so the validation section
 * (js/sections/validation.js) can render any type without knowing its internals:
 *   VALIDATION_CHART_TYPES[name] = {
 *     config: vc => Chart.js config object,
 *     legendHTML: vc => HTML string rendered above the canvas ('' to use the Chart.js legend),
 *   }
 * This file must load before the chart type files.
 */
const VALIDATION_TEXT_COLOR = '#333F52';
const VALIDATION_AXIS_TITLE_COLOR = '#074770';
const VALIDATION_GRID_COLOR = 'rgba(0,0,0,0.06)';

const VALIDATION_CHART_TYPES = {};

function validationAxisTitle(text, padding) {
  return { display: true, text, font: { size: 13, weight: '600' }, color: VALIDATION_AXIS_TITLE_COLOR, padding };
}

function hexToRgba(hex, alpha) {
  const n = parseInt(hex.replace('#', ''), 16);
  return `rgba(${(n >> 16) & 255}, ${(n >> 8) & 255}, ${n & 255}, ${alpha})`;
}
