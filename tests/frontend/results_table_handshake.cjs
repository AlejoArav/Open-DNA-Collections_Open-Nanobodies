// Regression: iframe loads before the host registers its message listener.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const script = fs.readFileSync(process.argv[2], 'utf8').match(/<script>([\s\S]*?)<\/script>/)[1];
const listeners = new Map(), timers = new Map(), messages = [];
let registered = false, timerId = 0, clickId = 0;
function element() {
  return {children: [], style: {}, scrollHeight: 50,
    replaceChildren() { this.children = []; },
    append(child) { this.children.push(child); },
    setAttribute() {}, dataset: {}};
}
const head = element(), body = element(), scroll = element();
const parent = {postMessage(message) { if (registered) messages.push(message); }};
const window = {parent, addEventListener(name, callback) { listeners.set(name, callback); }};
const document = {body: element(), documentElement: {style: {setProperty() {}}},
  createElement: element,
  querySelector(name) { return {thead: head, tbody: body, '.scroll': scroll}[name]; }};
vm.runInNewContext(script, {window, document,
  setInterval(callback) { timers.set(++timerId, callback); return timerId; },
  clearInterval(id) { timers.delete(id); },
  crypto: {randomUUID() { return `click-${++clickId}`; }}});
assert.equal(messages.length, 0, 'early readiness was dropped by the unregistered host');
assert.equal(timers.size, 1, 'readiness should retry while unacknowledged');
registered = true;
for (const callback of timers.values()) callback();
assert.equal(messages[0].type, 'streamlit:componentReady');
listeners.get('message')({source: parent, data: {type: 'streamlit:render',
  args: {columns: ['Name'], rows: [{part_key: 'BBF10K_000001', cells: ['Part one']}]}}});
assert.equal(timers.size, 0, 'the first render must stop readiness retries');
assert.equal(body.children.length, 1, 'the table should render after delayed registration');
body.children[0].onclick();
const clicked = messages.find(m => m.type === 'streamlit:setComponentValue');
assert.equal(clicked.value.part_key, 'BBF10K_000001');
assert.equal(clicked.value.event_id, 'click-1');
assert.equal(clicked.dataType, 'json');
console.log('Delayed host registration recovers and whole-row click succeeds');
