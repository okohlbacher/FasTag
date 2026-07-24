import type { ParamSpec } from '../../common/paramLayout'

export type ParamValue = string | boolean | string[]

interface Props {
  spec: ParamSpec
  value: ParamValue
  onChange: (v: ParamValue) => void
  onPickFile?: (spec: ParamSpec) => void
}

// One widget per declared type. The tool's restrictions drive the control:
// a string with `choices` becomes a select, a number carries its own min/max,
// so the UI cannot offer a value the CLI would reject.
export default function ParamField({ spec, value, onChange, onPickFile }: Props): JSX.Element {
  const label = (
    <label htmlFor={`p-${spec.name}`} title={spec.description}>
      {spec.name.replace(/_/g, ' ')}
    </label>
  )

  if (spec.type === 'bool') {
    return (
      <div className="field check">
        <input
          id={`p-${spec.name}`}
          type="checkbox"
          checked={value === true}
          onChange={(e) => onChange(e.target.checked)}
        />
        {label}
        <p className="help">{spec.description}</p>
      </div>
    )
  }

  if (spec.type === 'string-list') {
    const items = Array.isArray(value) ? value : []
    return (
      <div className="field">
        {label}
        <textarea
          id={`p-${spec.name}`}
          rows={Math.max(2, items.length + 1)}
          value={items.join('\n')}
          placeholder="one per line"
          onChange={(e) => onChange(e.target.value.split('\n'))}
        />
        <p className="help">{spec.description}</p>
      </div>
    )
  }

  if (spec.type === 'input-file' || spec.type === 'output-file') {
    return (
      <div className="field">
        {label}
        <div className="row tight">
          <input
            id={`p-${spec.name}`}
            type="text"
            value={String(value ?? '')}
            placeholder="(not set)"
            onChange={(e) => onChange(e.target.value)}
          />
          <button type="button" className="secondary slim" onClick={() => onPickFile?.(spec)}>
            Browse…
          </button>
        </div>
        <p className="help">{spec.description}</p>
      </div>
    )
  }

  if (spec.choices?.length) {
    return (
      <div className="field">
        {label}
        <select id={`p-${spec.name}`} value={String(value ?? '')} onChange={(e) => onChange(e.target.value)}>
          {spec.choices.map((c) => (
            <option key={c} value={c}>
              {c}
            </option>
          ))}
        </select>
        <p className="help">{spec.description}</p>
      </div>
    )
  }

  const numeric = spec.type === 'int' || spec.type === 'double'
  return (
    <div className="field">
      {label}
      <input
        id={`p-${spec.name}`}
        type={numeric ? 'number' : 'text'}
        step={spec.type === 'double' ? 'any' : 1}
        min={spec.min}
        max={spec.max}
        value={String(value ?? '')}
        onChange={(e) => onChange(e.target.value)}
      />
      <p className="help">{spec.description}</p>
    </div>
  )
}
