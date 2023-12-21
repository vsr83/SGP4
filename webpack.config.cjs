const path = require('path');

module.exports = {
  entry: './src/index.js',
  mode : 'production',
  output: {
    path: path.resolve(__dirname, 'dist'),
    filename: 'sgp4.js',
    library: {
        name: "sgp4",
        type: "umd",
    },
  },
};