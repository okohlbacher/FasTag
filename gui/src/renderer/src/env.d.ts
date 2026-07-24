/// <reference types="vite/client" />
import type { FastagApi } from '../../preload/index.d'

declare global {
  interface Window {
    fastag: FastagApi
  }
}
