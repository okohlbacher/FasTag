import React from 'react'
import { createRoot } from 'react-dom/client'
import './api' // installs window.fastag before App mounts
import App from './App'
import './index.css'

createRoot(document.getElementById('root')!).render(
  <React.StrictMode>
    <App />
  </React.StrictMode>
)
