const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const html = fs.readFileSync(process.argv[2], 'utf8');
const script = html.match(/<script>([\s\S]*?)<\/script>/)[1];
const listeners = new Map(), timers = new Map(), frames = [], messages = [], edits = [];
let registered = false, closed = 0, mounts = 0;
let width=0;
const root = {replaceChildren() {}, className: '', getBoundingClientRect() {return {width};}}, status = {textContent: ''};
const parent = {postMessage(message) { if (registered) messages.push(message); }};
const window = {addEventListener(name, fn) { listeners.set(name, fn); },
  ove: {createVectorEditor(node, options) {
    assert.ok(messages.some(m => m.type === 'streamlit:setFrameHeight'));
    assert.equal(options.readOnly, true);
    assert.equal(options.hideSingleImport, true);
    mounts++;
    return {updateEditor(data) { edits.push(data); }, close() {closed++;}};
  }}};
const document = {getElementById(id) {return id === 'viewer-root' ? root : status;},
  createElement() {return {style: {}};}};
vm.runInNewContext(script, {window, parent, document,
  setInterval(fn) {timers.set(1, fn); return 1;}, clearInterval(id) {timers.delete(id);},
  requestAnimationFrame(fn) {frames.push(fn);}});
timers.get(1)();
assert.equal(messages.length, 0, 'early ready dropped');
registered = true; timers.get(1)();
assert.equal(messages[0].type, 'streamlit:componentReady');
function render(hash, circular) {listeners.get('message')({source: parent,
  data: {type: 'streamlit:render', theme: {base:'dark'}, args: {payload: {content_hash: hash, sequenceData: {sequence:'ACGT',circular}}}}});}
render('linear', false);
assert.equal(timers.size, 0);
assert.equal(frames.length, 0, 'collapsed expander must not mount zero-width panels');
width=800;
listeners.get('load')();
assert.equal(mounts, 0, 'must wait for layout');
frames.shift()(); assert.equal(mounts, 0, 'one layout frame is insufficient');
frames.shift()(); assert.equal(mounts, 1);
assert.equal(root.className, '', 'viewer keeps its light theme inside a dark app');
assert.deepEqual(Array.from(edits[0].panelsShown[0], p => p.id), ['rail']);
render('linear', false); assert.equal(frames.length, 0, 'unchanged content avoids remount');
render('circle', true); frames.shift()(); frames.shift()();
assert.equal(closed, 1); assert.equal(mounts, 2);
assert.equal(edits[1].panelsShown[0][0].id, 'circular');
console.log('Viewer retries readiness, measures after two frames, stays read-only, and remounts on revision');
